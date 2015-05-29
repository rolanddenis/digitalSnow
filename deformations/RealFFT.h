#pragma once

#include <cstddef>
#include <stdexcept>
#include <new>

#include <complex> // To be included before fftw: see http://www.fftw.org/doc/Complex-numbers.html#Complex-numbers
#include <fftw3.h>

#include <DGtal/kernel/domains/HyperRectDomain.h>

#include <algorithm>

namespace DGtal
{

namespace
{

template <typename TFFTW>
struct FFTWComplexCast
  {
    static inline
    typename TFFTW::complex* apply( typename TFFTW::complex* ptr ) { return ptr; }

    static inline
    typename TFFTW::complex* apply( std::complex<typename TFFTW::real>* ptr ) { return reinterpret_cast<typename TFFTW::complex*>(ptr); }
  };


template <typename Real = double>
struct FFTWWrapper;

template <>
struct FFTWWrapper<double>
  {
    using size_t  = std::size_t;
    using real    = double;
    using complex = fftw_complex;
    using plan    = fftw_plan;
    using self    = FFTWWrapper<real>;

    static inline void*   malloc( size_t n )      { return fftw_malloc(n); }
    static inline void    free( void* p )         { fftw_free(p); }
    static inline void    execute( const plan p ) { fftw_execute(p); }
    static inline void    destroy_plan( plan p )  { fftw_destroy_plan(p); }


    template < typename C >
    static inline
    plan plan_dft_r2c( int rank, const int* n, real* in, C* out, unsigned flags )
      {
        return fftw_plan_dft_r2c(rank, n, in, FFTWComplexCast<self>::apply(out), flags);
      }

    template < typename C >
    static inline
    plan plan_dft_c2r( int rank, const int* n, C* in, real* out, unsigned flags )
      {
        return fftw_plan_dft_c2r(rank, n, FFTWComplexCast<self>::apply(in), out, flags);
      }

    template < typename C >
    static inline
    void execute_dft_r2c( const plan p, real* in, C* out )
      {
        fftw_execute_dft_r2c(p, in, FFTWComplexCast<self>::apply(out));
      }

    template < typename C >
    static inline
    void execute_dft_c2r( const plan p, C* in, real* out )
      {
        fftw_execute_dft_c2r(p, FFTWComplexCast<self>::apply(in), out);
      }
  };

} // anonymous namespace

/** Generic real-complex backward&forward FFT class
 * @tparam  TDomain Type of the domain over which the FFT will be performed.
 * @tparam  T       Values type.
 */
template <
  class TDomain, 
  typename T = double
>
class RealFFT;

/** Specialization for FFT over HyperRectDomain.
 * @tparam  TSpace  Type of the space.
 * @tparam  T       Values type.
 */
template <typename TSpace, typename T>
class RealFFT< HyperRectDomain<TSpace>, T >
  {
  private:
    using FFTW = FFTWWrapper<T>;

  public:
    using Space   = TSpace;
    using Domain  = HyperRectDomain<Space>;
    using Point   = typename Domain::Point;
    using Dimension = typename Domain::Dimension;
    using Real = T;
    using Complex = std::complex<Real>;
    static const Dimension dimension = Domain::dimension;

    /** Constructor.
     * @param aDomain The domain over which the transform will be performed.
     */
    RealFFT( Domain const& aDomain ) noexcept
        : mySpatialDomain{ aDomain }
        , mySpatialExtent{ mySpatialDomain.upperBound() - mySpatialDomain.lowerBound() + Point::diagonal(1) }
        , myFreqExtent{ mySpatialExtent / (Point::diagonal(1) + Point::base(0)) + Point::base(0) }
        , myFreqDomain{ Point::diagonal(0), myFreqExtent - Point::diagonal(1) }
        , myStorage( FFTW::malloc( sizeof(Complex) * myFreqDomain.size() ) )
      {}

    /// Destructor
    ~RealFFT()
      {
        FFTW::free( myStorage );
      }

    /** Checks if storage is valid.
     * @return true if there is an allocated storage, false otherwise.
     */
    bool isValid() const noexcept
      {
        return myStorage != nullptr;
      }

    /** Padding when using real datas. 
     *
     * @return the number of real values used as padding along the last dimension.
     *
     *  \see http://www.fftw.org/doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data 
     */
    inline  
    size_t getPadding() const noexcept
      {
        return 2*myFreqExtent[0] - mySpatialExtent[0];
      }

    /** Get mutable spatial storage.
     * \warning There is a padding at the end of the first dimension. \see getPadding
     */
    inline
    Real* getSpatialStorage() noexcept
      {
        return reinterpret_cast<Real*>(myStorage);
      }

    /** Get non-mutable spatial storage.
     * \warning There is a padding at the end of the first dimension. \see getPadding
     */
    inline
    Real const* getSpatialStorage() const noexcept
      {
        return reinterpret_cast<Real const*>(myStorage);
      }

    /// Get mutable frequential storage.
    inline
    Complex* getFreqStorage() noexcept
      {
        return reinterpret_cast<Complex*>(myStorage);
      }

    /// Get non-mutable frequential storage.
    inline
    Complex const* getFreqStorage() const noexcept
      {
        return reinterpret_cast<Complex const*>(myStorage);
      }
    
    /// Get spatial domain.
    inline Domain const& getSpatialDomain() const noexcept { return mySpatialDomain; }

    /// Get frequential domain.
    inline Domain const& getFreqDomain()    const noexcept { return myFreqDomain; }

    /// Get spatial domain extent.
    inline Point  const& getSpatialExtent() const noexcept { return mySpatialExtent; }
        
    /// Get frequential domain extent.
    inline Point  const& getFreqExtent()    const noexcept { return myFreqExtent; }

    /** Forward transformation (spatial -> frequential)
     *
     * @param flags Planner flags. \see http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags 
     */
    void forwardFFT( unsigned flags = FFTW_ESTIMATE )
      {
        typename FFTW::plan p;

        // Transform dimensions
        int n[dimension];
        for (size_t i = 0; i < dimension; ++i)
          n[dimension-i-1] = mySpatialExtent[i];

        // Create the plan for this transformation
        // Only FFTW_ESTIMATE flag preserve input.
        if ( flags & FFTW_ESTIMATE )
          {
            p = FFTW::plan_dft_r2c( dimension, n, getSpatialStorage(), getFreqStorage(), FFTW_ESTIMATE );
          }
        else
          {
            // Strategy to preserve input datas while creating DFT plan:
            // - Firstly, checks if a plan already exists for this dimensions.
            p = FFTW::plan_dft_r2c( dimension, n, getSpatialStorage(), getFreqStorage(), FFTW_WISDOM_ONLY | flags );

            // - Otherwise, create fake input to create the plan.
            if ( p == NULL )
              {
                void* tmp = FFTW::malloc( sizeof(Complex) * myFreqDomain.size() );
                if ( tmp == nullptr )  throw std::bad_alloc{};
                p = FFTW::plan_dft_r2c( dimension, n, reinterpret_cast<Real*>(tmp), reinterpret_cast<Complex*>(tmp), flags );
                FFTW::free(tmp);
              }
          }

        // We must have a valid plan now ...
        if ( p == NULL ) throw std::runtime_error("No valid DFT plan founded.");

        // Gogogo !
        FFTW::execute_dft_r2c( p, getSpatialStorage(), getFreqStorage() );

        // Destroying plan
        FFTW::destroy_plan( p );
      }

    /** Backward transformation (frequential -> spatial)
     *
     * @param flags Planner flags. \see http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags 
     */
    void backwardFFT( unsigned flags = FFTW_ESTIMATE )
      {
        typename FFTW::plan p;

        // Transform dimensions
        int n[dimension];
        for (size_t i = 0; i < dimension; ++i)
          n[dimension-i-1] = mySpatialExtent[i];

        // Create the plan for this transformation
        if ( flags & FFTW_ESTIMATE )
          {
            p = FFTW::plan_dft_c2r( dimension, n, getFreqStorage(), getSpatialStorage(), FFTW_ESTIMATE );
          }
        else
          {
            // Strategy to preserve input datas while creating DFT plan:
            // - Firstly, checks if a plan already exists for this dimensions.
            p = FFTW::plan_dft_c2r( dimension, n, getFreqStorage(), getSpatialStorage(), FFTW_WISDOM_ONLY | flags );
            
            // - Otherwise, create fake input to create the plan.
            if ( p == NULL )
              {
                void* tmp = FFTW::malloc( sizeof(Complex) * myFreqDomain.size() );
                if ( tmp == nullptr )  throw std::bad_alloc{};
                p = FFTW::plan_dft_c2r( dimension, n, reinterpret_cast<Complex*>(tmp), reinterpret_cast<Real*>(tmp), flags );
                FFTW::free(tmp);
              }
          }

        // We must have a valid plan now ...
        if ( p == NULL ) throw std::runtime_error("No valid DFT plan founded.");

        // Fire in the hole !
        FFTW::execute_dft_c2r( p, getFreqStorage(), getSpatialStorage() );

        // Destroying plan
        FFTW::destroy_plan( p );
      }

  private:
    const Domain  mySpatialDomain;  ///< Spatial domain (real).
    const Point   mySpatialExtent;  ///< Extent of the spatial domain.
    const Point   myFreqExtent;     ///< Extent of the frequential domain.
    const Domain  myFreqDomain;     ///< Frequential domain (complex).
          void*   myStorage;        ///< Storage.
    
  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
