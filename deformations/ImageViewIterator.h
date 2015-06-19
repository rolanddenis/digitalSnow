#pragma once

#include <type_traits>
#include <boost/iterator/iterator_facade.hpp>
#include "Linearizer.h"


namespace DGtal
{

template <
  typename TIterableClass
>
class ImageViewIterator
  : public boost::iterator_facade<
      ImageViewIterator<TIterableClass>,
      typename TIterableClass::Value,
      std::random_access_iterator_tag,
      decltype( ((TIterableClass*)nullptr)->dereference( TIterableClass::Point::diagonal(0), typename TIterableClass::Point::Coordinate(0) ) )
    >
  {
  public:
    using Self = ImageViewIterator<TIterableClass>;
    using IterableClass = TIterableClass;
    using Domain = typename IterableClass::Domain;
    using Point = typename Domain::Point;
    using Linearizer = DGtal::Linearizer<Domain, ColMajorStorage>; // Fixed but must be later set as template.
    using Reference = decltype( ((IterableClass*)nullptr)->dereference( Point::diagonal(0), typename Point::Coordinate(0) ) );

    /// Default constructor.
  ImageViewIterator()  
    : myIterableClassPtr(nullptr)
    {}

  /// Iterator from a point.
  ImageViewIterator( IterableClass* anIterableClassPtr, Domain const& aFullDomain, Domain const& aViewDomain, Point const& aPoint ) 
    : myIterableClassPtr( anIterableClassPtr )
    , myFullDomain{ aFullDomain }
    , myViewDomain{ aViewDomain }
    , myFullExtent( myFullDomain.upperBound() - myFullDomain.lowerBound() + Point::diagonal(1) )
    , myViewExtent( myViewDomain.upperBound() - myViewDomain.lowerBound() + Point::diagonal(1) )
    , myPoint{ aPoint }
    , myFullIndex( Self::Linearizer::getIndex( myPoint - myFullDomain.lowerBound(), myFullExtent ) )
    {}

  /// Copy constructor.
  template < typename TOtherIterableClass >
  ImageViewIterator( 
      ImageViewIterator<TOtherIterableClass> const& other,
      typename std::enable_if< std::is_convertible<TOtherIterableClass*, IterableClass*>::value >::type* = 0 
  )
    : myIterableClassPtr( other.myIterableClassPtr )
    , myFullDomain{ other.myFullDomain }
    , myViewDomain{ other.myViewDomain }
    , myFullExtent{ other.myFullExtent }
    , myViewExtent{ other.myViewExtent }
    , myPoint{ other.myPoint }
    , myFullIndex{ other.myFullIndex }
    {}

  /// Move constructor.
  template < typename TOtherIterableClass >
  ImageViewIterator( 
      ImageViewIterator<TOtherIterableClass> && other,
      typename std::enable_if< std::is_convertible<TOtherIterableClass*, IterableClass*>::value >::type* = 0 
  )
    : myIterableClassPtr( std::move(other.myIterableClassPtr) )
    , myFullDomain{ std::move(other.myFullDomain) }
    , myViewDomain{ std::move(other.myViewDomain) }
    , myFullExtent{ std::move(other.myFullExtent) }
    , myViewExtent{ std::move(other.myViewExtent) }
    , myPoint{ std::move(other.myPoint) }
    , myFullIndex{ std::move(other.myFullIndex) }
    {}

  /// Destructor.
  ~ImageViewIterator() 
    {}

  /// Copy assignment constructor.
  template < typename TOtherIterableClass >
  typename std::enable_if< std::is_convertible<TOtherIterableClass*, IterableClass*>::value, Self& >::type
  operator=( 
      ImageViewIterator<TOtherIterableClass> const& other
  ) 
    {
      myIterableClassPtr = other.myIterableClassPtr;
      myFullDomain = other.myFullDomain;
      myViewDomain = other.myViewDomain;
      myFullExtent = other.myFullExtent;
      myViewExtent = other.myViewExtent;
      myPoint = other.myPoint;
      myFullIndex = other.myFullIndex;
      return *this;
    }

  /// Move assignment constructor.
  template < typename TOtherIterableClass >
  typename std::enable_if< std::is_convertible<TOtherIterableClass*, IterableClass*>::value, Self& >::type
  operator=( 
      ImageViewIterator<TOtherIterableClass> && other
  ) 
    {
      myIterableClassPtr = std::move(other.myIterableClassPtr);
      myFullDomain = std::move(other.myFullDomain);
      myViewDomain = std::move(other.myViewDomain);
      myFullExtent = std::move(other.myFullExtent);
      myViewExtent = std::move(other.myViewExtent);
      myPoint = std::move(other.myPoint);
      myFullIndex = std::move(other.myFullIndex);
      return *this;
    }


public:


  /// begin iterator.
  explicit ImageViewIterator( IterableClass* anIterableClassPtr, Domain const& aFullDomain, Domain const& aViewDomain ) 
    : ImageViewIterator{ anIterableClassPtr, aFullDomain, aViewDomain, aViewDomain.lowerBound() }
    {}

  /// end iterator.
  explicit ImageViewIterator( IterableClass* anIterableClassPtr, Domain const& aFullDomain, Domain const& aViewDomain, bool /* last */ ) 
    : ImageViewIterator{ anIterableClassPtr, aFullDomain, aViewDomain, aViewDomain.upperBound() }
    {
      increment();
    }

  /// Return the point behind this iterator
  inline
  Point const& getPoint() const 
    {
      return myPoint;
    }

  /// Distance to a point
  inline
  std::ptrdiff_t distance_to( Point const& aPoint ) const 
    {
      return static_cast<std::ptrdiff_t>( Self::Linearizer::getIndex( aPoint - myPoint, myViewExtent ) );
    }

private:

  /// Friendship of interoperability.
  template <class> friend class ImageViewIterator;
  friend class boost::iterator_core_access;
  //friend class DistanceFunctor;
  

  /// Increment of one step.
  void increment() 
    {
      ++myFullIndex;
      ++myPoint[0];
      for ( std::size_t i = 1; i < Domain::dimension && myPoint[i-1] > myViewDomain.upperBound()[i-1]; ++i )
        {
          myPoint[i-1] = myViewDomain.lowerBound()[i-1];
          ++myPoint[i];
          std::size_t cum = myFullExtent[i-1] - myViewExtent[i-1];
          for ( std::size_t j = 0; j < i-1; ++j )
            cum *= myFullExtent[j];

          myFullIndex += cum;
        }
    }

  /// Decrement of one step.
  void decrement() 
    {
      --myFullIndex;
      --myPoint[0];
      for ( std::size_t i = 1; i < Domain::dimension && myPoint[i-1] < myViewDomain.lowerBound()[i-1]; ++i )
        {
          myPoint[i-1] = myViewDomain.upperBound()[i-1];
          --myPoint[i];
          std::size_t cum = myFullExtent[i-1] - myViewExtent[i-1];
          for ( std::size_t j = 0; j < i-1; ++j )
            cum *= myFullExtent[j];

          myFullIndex -= cum;
        }
    }

  /// Equality.
  inline
  bool equal( Self const& other ) const 
    {
      return myFullIndex == other.myFullIndex;
    }

  /// Constant dereference.
    inline
    Reference dereference() const
      {
        return myIterableClassPtr->dereference( myPoint, myFullIndex );
      }

    /// Distance to other iterator.
    inline
    std::ptrdiff_t distance_to( Self const& other ) const 
      {
        return distance_to( other.myPoint );
      }

    /// Advance by n steps. Not very efficient implementation ...
    void advance( std::ptrdiff_t n ) 
      {
        const auto pos = Self::Linearizer::getIndex( myPoint, myViewDomain.lowerBound(), myViewExtent );
        myPoint = Self::Linearizer::getPoint( pos + n, myViewDomain.lowerBound(), myViewExtent );
        myFullIndex = Self::Linearizer::getIndex( myPoint, myFullDomain.lowerBound(), myFullExtent ); 
      }

  private:
    IterableClass* myIterableClassPtr;
    Domain myFullDomain;
    Domain myViewDomain;
    Point myFullExtent;
    Point myViewExtent;
    Point myPoint;
    typename Point::Coordinate myFullIndex;
    
  };




} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
