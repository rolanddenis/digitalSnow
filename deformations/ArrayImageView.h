#pragma once

#include <boost/assert.hpp>
#include <boost/concept/assert.hpp>
#include <boost/iterator/iterator_concepts.hpp>
#include <iterator>

#include "ImageViewIterator.h"
#include "IteratorFacade.h"
#include "Linearizer.h"

namespace DGtal
{

/** Image view for c-style array.
 *
 * This creates an image (concepts::CImage compatible) given a pointer to his
 * allocated memory and two domains:
 * - the definition (full) domain whose size is equal to the allocated size.
 * - the viewable domain, a subset of the full-domain, on which the image is accessible.
 *
 * The available iterators can return the corresponding point and are faster than using
 * an iterator over the domain (see ImageViewIterator). Reverse iterators and ranges are 
 * defined in the inherited class IteratorFacade.
 *
 * @warning The allocated memory must be column-major ordered.
 * @warning The domain is supposed to be an HyperRectDomain.
 *
 * @tparam TDomain  Type of the domain (must be an HyperRectDomain).
 * @tparam TValue   Type of the stored values.
 */
template < 
  typename TDomain,
  typename TArrayIterator
>
class ArrayImageView
    : public IteratorFacade< ArrayImageView<TDomain,TArrayIterator> >
  {

  // Check Random-access iterator concept on TIterator
  BOOST_CONCEPT_ASSERT( (boost_concepts::RandomAccessTraversalConcept<TArrayIterator>) );

  public:
    // Typedefs
    using Self = ArrayImageView<TDomain, TArrayIterator>;
    using ArrayIterator = TArrayIterator;
    using Value = typename std::iterator_traits<ArrayIterator>::value_type;
    using Reference = typename std::iterator_traits<ArrayIterator>::reference;
    using ConstReference = const Reference;

    using Domain = TDomain;
    using Point = typename Domain::Point;
    using Linearizer = DGtal::Linearizer<Domain, ColMajorStorage>; ///< Linearization of the points.

    // Iterators & Ranges
    template <class> friend class ImageViewIterator;
    using Iterator = typename IteratorTraits<Self>::Iterator; ///< Mutable iterator.
    using ConstIterator = typename IteratorTraits<Self>::ConstIterator; ///< Constant iterator.

    /** Default constructor.
     *
     * Empty allocated memory on empty domains.
     */
    ArrayImageView()
      : myArrayIterator{nullptr}
      , myFullDomain{}
      , myViewDomain{}
      {}

    /// Constructor from storage, full domain and viewable domain.
    ArrayImageView( ArrayIterator anArrayIterator, Domain const& aFullDomain, Domain const& aViewDomain )
        : myArrayIterator(anArrayIterator)
        , myFullDomain{ aFullDomain }
        , myViewDomain{ aViewDomain }
      {
        BOOST_ASSERT_MSG(
               aFullDomain.lowerBound().isLower( aViewDomain.lowerBound() )
            && aFullDomain.upperBound().isUpper( aViewDomain.upperBound() ),
            "The viewable domain must be included into the full domain."
        );
      }

    /** Constructor from storage and full domain.
     *
     * The viewable domain is then the full domain.
     */
    ArrayImageView( ArrayIterator anArrayIterator, Domain const& aFullDomain )
        : ArrayImageView( anArrayIterator, aFullDomain, aFullDomain )
      {
      }

    /// Copy constructor with other viewable domain.
    ArrayImageView( Self const& other, Domain const& aViewDomain )
        : ArrayImageView( other.myArrayIterator, other.myFullDomain, aViewDomain )
      {}
   
    /**
     * @return the image viewable domain.
     */
    inline
    Domain domain() const
      {
        return myViewDomain;
      }

    /**
     * @return the full domain where the allocated memory is defined.
     */
    inline
    Domain fullDomain() const
      {
        return myFullDomain;
      }
    
    /**
     * @return the value given a point lying inside the full domain.
     */
    inline
    Value getValue( Point const& aPoint ) const
      {
        BOOST_ASSERT_MSG(
            myFullDomain.isInside(aPoint),
            "The point is outside the full domain."
        );

        return myArrayIterator[ Linearizer::getIndex(aPoint, myFullDomain) ];
      }

    /// Set a value given a point lying inside the full domain.
    inline
    void setValue( Point const& aPoint, Value aValue )
      {
        BOOST_ASSERT_MSG(
            myFullDomain.isInside(aPoint),
            "The point is outside the full domain."
        );

        myArrayIterator[ Linearizer::getIndex(aPoint, myFullDomain) ] = aValue;
      }

    /**
     * @return the value given a point lying inside the full domain.
     */
    inline
    Value operator() ( Point const& aPoint ) const
      {
        return getValue(aPoint);
      }

    /**
     * @return an mutable iterator pointing to the lower bound of the viewable domain.
     */
    inline
    Iterator begin()
      {
        return Iterator{ this, myFullDomain, myViewDomain };
      }
    
    /**
     * @return a constant iterator pointing to the lower bound of the viewable domain.
     */
    inline
    ConstIterator begin() const
      {
        return ConstIterator{ this, myFullDomain, myViewDomain };
      }
    
    /**
     * @return a constant iterator pointing to the lower bound of the viewable domain (C++11).
     */
    inline
    ConstIterator cbegin() const
      {
        return ConstIterator{ this, myFullDomain, myViewDomain };
      }

    /**
     * @return an mutable iterator pointing after the upper bound of the viewable domain.
     */
    inline
    Iterator end()
      {
        return Iterator{ this, myFullDomain, myViewDomain, true };
      }
    
    /**
     * @return a constant iterator pointing after the upper bound of the viewable domain.
     */
    inline
    ConstIterator end() const
      {
        return ConstIterator{ this, myFullDomain, myViewDomain, true };
      }
    
    /**
     * @return a constant iterator pointing after the upper bound of the viewable domain (C++11).
     */
    inline
    ConstIterator cend() const
      {
        return ConstIterator{ this, myFullDomain, myViewDomain, true };
      }
      

  public: // Should be private since ImageViewIterator is a friend but g++ 4.9.1 don't care ... (no prob with clang++)

    /// Dereference of a mutable iterator.
    inline
    Reference dereference( Point const& /* aPoint */, typename Point::Coordinate aFullIndex )
      {
        BOOST_ASSERT_MSG(
            aFullIndex >= 0 && static_cast<typename Domain::Size>(aFullIndex) < myFullDomain.size(),
            "linearized index out of bounds !"
        );
        return myArrayIterator[aFullIndex];
      }

    /// Dereference of a constant iterator.
    inline
    ConstReference dereference( Point const& /* aPoint */, typename Point::Coordinate aFullIndex ) const
      {
        BOOST_ASSERT_MSG(
            aFullIndex >= 0 && static_cast<typename Domain::Size>(aFullIndex) < myFullDomain.size(),
            "linearized index out of bounds !"
        );
        return myArrayIterator[aFullIndex];
      }

  private:
    ArrayIterator myArrayIterator; ///< Pointer to the allocated memory.
    Domain myFullDomain;  ///< Definition (full) domain.
    Domain myViewDomain;  ///< Viewable domain.
  };

  /** Iterator traits specialized for CArrayImageView.
   *
   * \see ImageViewIterator
   */
  template <
    typename TDomain,
    typename TArrayIterator
  >
  class IteratorTraits< ArrayImageView<TDomain, TArrayIterator> >
    {
    public:
      using Self = ArrayImageView<TDomain, TArrayIterator>;
      using Iterator = ImageViewIterator<Self>; ///< Mutable iterator.
      using ConstIterator = ImageViewIterator<const Self>; ///< Constant iterator.

      /** Functor that return the distance between the domain lower bound and a given point.
       *
       * \see SimpleRandomAccessRangeFromPoint and SimpleRandomAccessConstRangeFromPoint.
       */
      class DistanceFunctor
        {
        public:
          using Domain = typename Self::Domain;
          using Point = typename Self::Point;
          using Difference = typename Self::Difference;

          DistanceFunctor( Self const* anImageView )
            : myDomain( anImageView->domain() )
            {}

          Difference operator() ( Point const& aPoint ) const
            {
              BOOST_ASSERT_MSG(
                  myDomain.isInside(aPoint),
                  "The point is outside the domain !"
              );
              return Linearizer<Domain, ColMajorStorage>::getIndex( aPoint, myDomain );
            }

        private:
          const Domain myDomain; ///< Stored domain to avoid iterator corruption if domain changed.
        };

    };

} // namespace DGtal


/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
