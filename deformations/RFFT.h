#pragma once

#include <cstddef>

#include <complex> // To be includes before fftw:  http://www.fftw.org/doc/Complex-numbers.html#Complex-numbers
#include <fftw3.h>

#include <DGtal/kernel/domains/HyperRectDomain.h>

#include <algorithm>

namespace DGtal
{

template <typename Real = double>
class FFTW;

template <>
class FFTW<double>
  {
    using size_t  = std::size_t;
    using real    = double;
    using complex = fftw_complex;
    using plan    = fftw_plan;
    static inline void*   malloc( size_t n )      const { fftw_malloc(n); }
    static inline void    free( void* p )         const { fftw_free(p); }
    static inline void    execute( const plan p ) const { fftw_execute(p); }
    static inline void    destroy_plan( plan p )  const { fftw_destroy_plan(p); }
    
    static inline
    plan plan_dft_r2c( int rank, const int* n, real* in, complex* out, unsigned flags ) const
      {
        return fftw_plan_dft_r2c(rank, n, in, out, flags);
      }

    static inline
    plan plan_dft_c2r( int rank, const int* n, complex* in, real* out, unsigned flags ) const
      {
        return fftw_plan_dft_c2r(rank, n, in, out, flags);
      }

    static inline
    void execute_dft_r2c( const plan p, real* in, complex* out ) const
      {
        return fftw_execute_dft_r2c(p, in, out);
      }

    static inline
    void execute_dft_c2r( const plan p, complex* in, real* out ) const
      {
        return fftw_execute_dft_c2r(p, in, out);
      }
  };


template <
  class TDomain, 
  typename T = double
>
class RealFFT;

template <typename TSpace>
class RealFFT< HyperRectDomain<TSpace> >
  {
  public:
    using Space = TSpace;
    using Domain = HyperRectDomain<Space>;
    using Point = Domain::Point;
    using Dimension = Domain::Dimension;
    using Real = T;
    using FFTW = FFTW<T>;
    using Complex = FFTW::complex;
    static const Dimension dimension = Domain::dimension;

    RealFFT( Domain const& aDomain )
        : myDomain(aDomain)
      {}

    

  private:
    Domain myDomain;
    
  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
