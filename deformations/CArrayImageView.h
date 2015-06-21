#pragma once

#include <boost/static_assert.hpp>

#include "ImageViewIterator.h"
#include "IteratorFacade.h"
#include "Linearizer.h"

namespace DGtal
{

template < 
  typename TDomain,
  typename TValue
>
class CArrayImageView
    : public IteratorFacade< CArrayImageView<TDomain,TValue> >
  {

  public:
    // Typedefs
    using Self = CArrayImageView<TDomain, TValue>;
    using Value = TValue;
    using Domain = TDomain;
    using Point = typename Domain::Point;
    using Linearizer = DGtal::Linearizer<Domain, ColMajorStorage>;

    // Iterators & Ranges
    template <class> friend class ImageViewIterator;
    using Iterator = typename IteratorTraits<Self>::Iterator;
    using ConstIterator = typename IteratorTraits<Self>::ConstIterator;

    /// Default constructor.
    CArrayImageView()
      : myStorage{nullptr}
      {}

    /// Constructor from storage, full domain and viewable domain.
    CArrayImageView( Value* aStorage, Domain const& aFullDomain, Domain const& aViewDomain )
        : myStorage(aStorage)
        , myFullDomain{ aFullDomain }
        , myViewDomain{ aViewDomain }
      {
        BOOST_ASSERT_MSG(
               aFullDomain.lowerBound().isLower( aViewDomain.lowerBound() )
            && aFullDomain.upperBound().isUpper( aViewDomain.upperBound() ),
            "The viewable domain must be included into the full domain."
        );
      }

    /// Constructor from storage and full domain.
    CArrayImageView( Value* aStorage, Domain const& aFullDomain )
        : CArrayImageView( aStorage, aFullDomain, aFullDomain )
      {
      }

    /// Copy constructor with other viewable domain.
    CArrayImageView( Self const& other, Domain const& aViewDomain )
        : CArrayImageView( other.myStorage, other.myFullDomain, aViewDomain )
      {}
   
    /// Get image domain.
    inline
    Domain domain() const
      {
        return myViewDomain;
      }
    
    /// Get a value
    inline
    Value getValue( Point const& aPoint ) const
      {
        return myStorage[ Linearizer::getIndex(aPoint, myFullDomain) ];
      }

    /// Set a value
    inline
    void setValue( Point const& aPoint, Value aValue )
      {
        myStorage[ Linearizer::getIndex(aPoint, myFullDomain) ] = aValue;
      }

    /// Get a value
    inline
    Value operator() ( Point const& aPoint ) const
      {
        return getValue(aPoint);
      }

    /// Iterators
    inline
    Iterator begin()
      {
        return Iterator{ this, myFullDomain, myViewDomain };
      }

    inline
    Iterator end()
      {
        return Iterator{ this, myFullDomain, myViewDomain, true };
      }

    inline
    ConstIterator cbegin() const
      {
        return ConstIterator{ this, myFullDomain, myViewDomain };
      }

    inline
    ConstIterator cend() const
      {
        return ConstIterator{ this, myFullDomain, myViewDomain, true };
      }
      

  private:

    inline
    Value& dereference( Point const& /* aPoint */, typename Point::Coordinate aFullIndex )
      {
        return myStorage[aFullIndex];
      }

    inline
    Value const& dereference( Point const& /* aPoint */, typename Point::Coordinate aFullIndex ) const
      {
        return myStorage[aFullIndex];
      }

  private:
    Value* myStorage;
    Domain myFullDomain;
    Domain myViewDomain;
  };

  /// Iterator traits
  template <
    typename TDomain,
    typename TValue
  >
  class IteratorTraits< CArrayImageView<TDomain, TValue> >
    {
    public:
      using Self = CArrayImageView<TDomain, TValue>;
      using Iterator = ImageViewIterator<Self>;
      using ConstIterator = ImageViewIterator<const Self>;

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
            return Linearizer<Domain, ColMajorStorage>::getIndex( aPoint, myDomain );
          }

      private:
        const Domain myDomain;
      };

    };

} // namespace DGtal


/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
