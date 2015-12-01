#pragma once

#include <type_traits>
#include <utility>
#include <iterator>
#include <functional>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/concept/assert.hpp>

#include <DGtal/kernel/domains/CDomain.h>

namespace DGtal
{

/* Transform a point-dependant functor into a constant image.
 *
 * @tparam TDomain  Domain type.
 * @tparam TValue   Value type returned by the functor.
 *
 * The functor must accept a point and return a value whose type is \a Value.
 * Since the functor will be stored by copy, prefer a lightweight type.
 */
template <
  typename TDomain,
  typename TValue
>
class FunctorConstImage
  {
    BOOST_CONCEPT_ASSERT(( DGtal::concepts::CDomain<TDomain> ));

public:
    // DGtal types
    using Self      = FunctorConstImage<TDomain, TValue>;
    using Domain    = TDomain;
    using Point     = typename Domain::Point;
    using Vector    = typename Domain::Vector;
    using Integer   = typename Domain::Integer;
    using Size      = typename Domain::Size;
    using Dimension = typename Domain::Dimension;
    using Vertex    = Point;
    using Value     = TValue;
    using Functor   = std::function<Value(Point)>;
    
    using ConstIterator = boost::transform_iterator< std::reference_wrapper<const Functor>, typename Domain::ConstIterator >;
    using ConstReverseIterator = std::reverse_iterator< ConstIterator >;
    class ConstRange;

    BOOST_STATIC_CONSTANT( Dimension, dimension = Domain::Space::dimension );

public:
    /** Constructor
     * @param aDomain   The domain of the image.
     * @param aFunctor  The functor taking point as parameter.
     */
    template < class TFunctor >
    FunctorConstImage( Domain const& aDomain, TFunctor && aFunctor )
      : myDomain( aDomain )
      , myFunctor( std::forward<TFunctor>(aFunctor) )
    {
    }

    /**
     * @return the associated domain.
     */
    inline
    Domain const& domain() const
      {
        return myDomain;
      }

    /** Gets the value of the functor for the given point.
     * @param aPoint the point.
     * @return the value at \a aPoint.
     */
    inline
    Value operator() ( Point const& aPoint ) const
      {
        return myFunctor( aPoint );
      }
    
    /**
     * @return a constant range over this image.
     */
    inline
    ConstRange constRange() const 
      {
        return ConstRange( *this );
      }

public:
    /// Constant range
    class ConstRange
      {
    public:
        ConstRange( Self const& aFunctorConstImage )
          : myFunctorConstImage( aFunctorConstImage )
        {}

        using ConstIterator = Self::ConstIterator;
        using ConstReverseIterator = Self::ConstReverseIterator;
        using Point = Self::Point;

        inline ConstIterator begin()  const { return { myFunctorConstImage.myDomain.begin(), myFunctorConstImage.myFunctor }; }
        inline ConstIterator begin( Point const& aPoint ) const { return { myFunctorConstImage.myDomain.begin(aPoint), myFunctorConstImage.myFunctor }; }
        inline ConstIterator end()    const { return { myFunctorConstImage.myDomain.end(), myFunctorConstImage.myFunctor }; }
        
        inline ConstReverseIterator rbegin()  const { return ConstReverseIterator( end() ); }
        inline ConstReverseIterator rbegin( Point const& aPoint ) const { return ConstReverseIterator( ++begin(aPoint) ); }
        inline ConstReverseIterator rend()    const { return ConstReverseIterator( begin() ); }

    private:
        Self const& myFunctorConstImage;
      };


private:
    Domain  myDomain;   ///< The image domain.
    Functor myFunctor;  ///< The functor that generates the image.
  };

/** FunctorConstImage construction helper.
 *
 * @tparam  TDomain   The domain type (auto-deduced).
 * @tparam  TFunctor  The functor type (auto-deduced).
 * @param   aDomain   The image domain.
 * @param   aFunctor  The functor that generates the image.
 * @return an instance of the appropriate FunctorConstImage type.
 */
template <
  typename TDomain,
  typename TFunctor
>
FunctorConstImage< TDomain, typename std::result_of<TFunctor(typename TDomain::Point)>::type >
makeFunctorConstImage( TDomain const& aDomain, TFunctor && aFunctor )
  {
    return { aDomain, std::forward<TFunctor>(aFunctor) };
  }


} // namespace DGtal
