#pragma once

#include <boost/concept_check.hpp>
#include <cstddef>
#include <vector>

#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/domains/CDomain.h>
#include "ValueApproximations.h"
#include "NoBoundingBox.h"
#include "ImageView.h"

namespace 
{

  /**
   * Templated static structure for linearization of the coordinates of a point.
   *
   * @tparam TDomain    Type of the domain.
   * @tparam dimension  Actual working dimension (recursive process).
   */
  template < typename TDomain, std::size_t dimension = TDomain::dimension >
  struct linearizer
    {
      using Point = typename TDomain::Point; ///< Type of a point for which we want the linearized index.
      using Size =  typename TDomain::Size;  ///< Type of an index.

      /**
       * Return the linearized index from the coordinates of a point.
       *
       * @param aPoint      The point.
       * @param lowerBound  The lowerBound of the domain.
       * @param extent      The extent of the domain.
       */
      static inline
      Size apply( Point const& aPoint, Point const& lowerBound, Point const& extent ) noexcept
        {
          return 
              ( aPoint[TDomain::dimension - dimension] - lowerBound[TDomain::dimension - dimension] ) 
            + extent[TDomain::dimension - dimension] * linearizer<TDomain, dimension-1>::apply( aPoint, lowerBound, extent );
        }
    };

  /**
   * Specialization of the structure linearizer for the one dimensional case.
   *
   * It is actually used as a terminate condition of the recursive process.
   */
  template < typename TDomain >
  struct linearizer< TDomain, 1 >
    {
      using Point = typename TDomain::Point;
      using Size  = typename TDomain::Size;

      static inline
      Size apply( Point const& aPoint, Point const& lowerBound, Point const& /* extent */ ) noexcept
        {
          return aPoint[TDomain::dimension-1] - lowerBound[TDomain::dimension-1];
        }
    };

  /**
   * Type traits to get the value type of a container.
   *
   * For most containers, the typedef value_type is what we want but others containers
   * like std::map of DGtal::LabelledMap exhibits instead the mapped_type typedef 
   * ( value_type represents, in this cases, the pair key|value )
   *
   * @tparam TContainer Type of the container
   */
  template < typename TContainer >
  struct ValueType
    {
      typedef typename TContainer::value_type type; ///< The type of the contained data.
    };

  /**
   * Specialization of ValueType for LabelledMap.
   */
  template < typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M >
  struct ValueType< DGtal::LabelledMap<TData, L, TWord, N, M> >
    {
      typedef TData type;
    };


} // namespace

namespace DGtal
{

/**
 * Multiple images container with approximation and bounding box capabilities.
 *
 * @tparam TDomain    The domain of the images.
 * @tparam TContainer The container used to stored, at each point of the domain, the values of the images.
 * @tparam TApproximation The predicate used to approximate the values.
 * @tparam TBoundingBox   The type of bounding box for the non-approximated values.
 */
template <
  typename TDomain,
  typename TContainer,
  typename TApproximation = approximations::NoValueApproximation< typename ValueType<TContainer>::type >,
  typename TBoundingBox = NoBoundingBox<TDomain>
>
class ApproximatedMultiImage;

/**
 * Specialization of ApproximatedMultiImage for LabelledMap container.
 */
template <
  typename TDomain,
  typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M,
  typename TApproximation,
  typename TBoundingBox
>
class ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox >
  {
  public:

    BOOST_CONCEPT_ASSERT( (concepts::CDomain<TDomain>) );
    BOOST_CONCEPT_ASSERT( (concepts::CValueApproximation<TApproximation>) );
    static_assert( std::is_same< typename TApproximation::value_type, TData>::value, "Container and approximation must have same value_type");
    /// \todo Concept check for bounding box

    // Typedefs
    template <typename> class ApproximatedReference;

    // DGtal typedefs
    using Self          = ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox >;
    using Domain        = TDomain; ///< Type of the domain.
    using Container     = DGtal::LabelledMap<TData, L, TWord, N, M>; ///< Type of the values container.
    using Approximation = TApproximation; ///< Type of the approximation.
    using BoundingBox   = TBoundingBox;   ///< Type of the bounding box.
    using Label         = typename Container::Label;  ///< Type of the label used to identify an image.
    using Value         = TData;  ///< Type of the stored values.
    using Point         = typename TDomain::Point;  ///< Type of a point in the domain.
    using Size          = typename TDomain::Size;   ///< Type of the linearized index of a point.
    using ValueConstIterator = typename Container::ConstIterator;   ///< Type of the const-iterator over the values associated to a given point.
    using FullImage     = ImageView< Self, image_view::FullDomain >; ///< Full image view.
    using ConstFullImage = ImageView< const Self, image_view::FullDomain >; ///< Full image constant view;
    using BBImage       = ImageView< Self, image_view::BoundingBoxAsDomain >; ///< Image view from bounding box.
    using ConstBBImage  = ImageView< const Self, image_view::BoundingBoxAsDomain >; ///< Image constant view from bounding box.


    // STL typedefs
    // ...
    
    /**
     * Constructor.
     */
    ApproximatedMultiImage( Domain const& aDomain, Approximation const& anApprox = Approximation{} )
      : myDomain{aDomain}, myImages{aDomain.size()}, 
        myApproximation{anApprox}, myBoundingBoxes{L, BoundingBox{aDomain}},
        myExtent{ aDomain.upperBound() - aDomain.lowerBound() + Point::diagonal(1) }
    {}

    /**
     * Get a value.
     *
     * @param aPoint The point of the domain.
     * @param aLabel The label of the image from which to get the value.
     */
    Value getValue( Point const& aPoint, Label aLabel ) const noexcept
      {
        Container const& values = myImages[ linearized(aPoint) ];
        
        if ( values.count(aLabel) > 0 )
          {
            return values.fastAt(aLabel);
          } 
        else
          {
            return myApproximation.default_value;
          }
      }
        
    /**
     * Set a value.
     *
     * @param aPoint  The point of the domain.
     * @param aLabel  The label of the image in which to set the value.
     * @param aValue  The value to be set.
     */
    void setValue( Point const& aPoint, Label aLabel, Value aValue ) noexcept
      {
        Container & values = myImages[ linearized(aPoint) ];

        if ( myApproximation.eval( aValue ) )
          {
            const auto it = values.find(aLabel);
            if ( it != values.end() )
              {
                values.erase(it);
                myBoundingBoxes[aLabel].removePoint(aPoint);
              }
          }
        else
          {
            if ( values.count(aLabel) > 0 )
              {
                values.fastAt(aLabel) = aValue;
              }
            else
              {
                values[aLabel] = aValue;
                myBoundingBoxes[aLabel].addPoint(aPoint);
              }
          }
      }

    /**
     * Return the domain of the images.
     */
    inline
    Domain const& domain() const noexcept
      {
        return myDomain;
      }
    
    /**
     * Return the bounding box (with buffer) of an image.
     * @param aLabel  Label of the image.
     * @param buffer  Buffer around the initial bounding box.
     */
    inline
    Domain getBoundingBox( Label aLabel, Point const& buffer = Point::diagonal(0) ) const noexcept
      {
        return myBoundingBoxes[ aLabel ].getBoundingBox( buffer );
      }
    
    /**
     * Return the value container associated to a point.
     * \todo mutable version.
     *
     * @param aPoint The point.
     */
    inline
    Container const& operator() ( Point const& aPoint ) const noexcept
      {
        return myImages[ linearized(aPoint) ];
      }

    /**
     * Return a full-domain view of the image with the given label
     * @param aLabel The image label.
     */
    inline
    FullImage operator[] ( Label aLabel ) noexcept
      {
        return { *this, aLabel };
      }

    /**
     * Return a full-domain constant view of the image with the given label
     * @param aLabel The image label.
     */
    inline
    ConstFullImage operator[] ( Label aLabel ) const noexcept
      {
        return { *this, aLabel };
      }

    /**
     * Return an image view restricted to his bounding box (with buffer).
     * @param aLabel  The image label.
     * @param buffer  The buffer around the bounding box.
     */
    inline
    BBImage getBBImage( Label aLabel, Point buffer = Point::diagonal(0) )
      {
        BBImage image{ *this, aLabel };
        image.buffer() = buffer;
        return image;
      }
    
    /**
     * Return an image constant-view restricted to his bounding box (with buffer).
     * @param aLabel  The image label.
     * @param buffer  The buffer around the bounding box.
     */
    inline
    ConstBBImage getBBImage( Label aLabel, Point buffer = Point::diagonal(0) ) const
      {
        ConstBBImage image{ *this, aLabel };
        image.buffer() = buffer;
        return image;
      }

  private:

    /**
     * Return the linearized index of a point.
     *
     * @param aPoint The point.
     */
    inline
    Size linearized( Point const& aPoint ) const
      {
        return linearizer<Domain>::apply(
            aPoint,
            myDomain.lowerBound(),
            myExtent
        );
      }


  private:
    Domain myDomain;
    std::vector<Container> myImages;
    Approximation myApproximation;
    std::vector<BoundingBox> myBoundingBoxes;
    Point myExtent;
    
  };

} // namespace DGtal


/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

