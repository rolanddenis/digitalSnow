#pragma once

#include <boost/concept_check.hpp>
#include <cstddef>
#include <vector>

#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/domains/CDomain.h>
#include "ValueApproximations.h"
#include "NoBoundingBox.h"

namespace 
{

  template < typename TDomain, std::size_t dimension = TDomain::dimension >
  struct linearizer
    {
      using Point = typename TDomain::Point;
      using Size =  typename TDomain::Size;

      static inline
      Size apply( Point const& aPoint, Point const& lowerBound, Point const& extent ) noexcept
        {
          return 
              ( aPoint[TDomain::dimension - dimension] - lowerBound[TDomain::dimension - dimension] ) 
            + extent[TDomain::dimension - dimension] * linearizer<TDomain, dimension-1>::apply( aPoint, lowerBound, extent );
        }
    };

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

  template < typename TContainer >
  struct ValueType
    {
      typedef typename TContainer::value_type type;
    };

  template < typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M >
  struct ValueType< DGtal::LabelledMap<TData, L, TWord, N, M> >
    {
      typedef TData type;
    };


} // namespace

namespace DGtal
{

template <
  typename TDomain,
  typename TContainer,
  typename TApproximation = approximations::NoValueApproximation< typename ValueType<TContainer>::type >,
  typename TBoundingBox = NoBoundingBox<TDomain>
>
class ApproximatedMultiImage;

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

    using Self          = ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox >;
    using Domain        = TDomain;
    using Container     = DGtal::LabelledMap<TData, L, TWord, N, M>;
    using Approximation = TApproximation;
    using BoundingBox   = TBoundingBox;
    using Label         = typename Container::Label;
    using Value         = TData;
    using Point         = typename TDomain::Point;
    using Size          = typename TDomain::Size;
    
    /**
     * Constructor.
     */
    ApproximatedMultiImage( Domain const& aDomain, Approximation const& anApprox = Approximation{} )
      : myDomain{aDomain}, myImages{aDomain.size()}, 
        myApproximation{anApprox}, myBoundingBoxes{L, BoundingBox{aDomain}},
        myExtent{ aDomain.upperBound() - aDomain.lowerBound() + Point::diagonal(1) }
    {}

    /**
     * Get a value
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
     * Set a value
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

  private:

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

