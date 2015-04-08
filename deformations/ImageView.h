#pragma once

#include <cstddef>
#include <array>
#include <type_traits>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include <DGtal/images/ImageContainerBySTLVector.h> // For conversion purpose

namespace DGtal
{

namespace image_view
{

/**
 * Policy that returns the full domain
 */
template < typename TImageView, typename TDomain >
struct FullDomain
  {
    inline
    TDomain domain() const
      {
        return static_cast<TImageView const*>(this)->myMultiImage->domain();
      }
  };

/**
 * Policy that returns the bounding box (with buffer) as the image domain.
 */
template < typename TImageView, typename TDomain >
class BoundingBoxAsDomain
  {
  public:
    using Point = typename TDomain::Point;

    inline
    TDomain domain() const
      {
        TImageView const& image = * static_cast<TImageView const*>(this);
        return image.myMultiImage->getBoundingBox( image.myLabel, myBuffer );
      }

    inline
    Point & buffer() 
      {
        return myBuffer;
      }

    inline
    Point const& buffer() const
      {
        return myBuffer;
      }

  private:
    Point myBuffer;
  };

} // namespace image_view


/**
 * Image View for ApproximatedMultiImage (but not only ?)
 */
template <
  typename TMultiImage,
  template<typename,typename> class TDomainPolicy = image_view::FullDomain
>
class ImageView
    : public TDomainPolicy< ImageView<TMultiImage, TDomainPolicy>, typename TMultiImage::Domain >
  {
  public:
    // Typedefs
    using Self        = ImageView<TMultiImage, TDomainPolicy>;
    using MultiImage  = TMultiImage;
    using Point       = typename MultiImage::Point;
    using Label       = typename MultiImage::Label;
    using Domain      = typename MultiImage::Domain;
    using Value       = typename MultiImage::Value;
    using Reference   = typename MultiImage::Reference;

    // Iterators & Ranges
    template <class> class MyIterator;
    class DistanceFunctor;
    using ConstIterator   = MyIterator<Self const>;
    using Iterator        = typename std::conditional< std::is_const<TMultiImage>::value, ConstIterator, MyIterator<Self> >::type;
    using ReverseIterator = boost::reverse_iterator<Iterator>;
    using ConstReverseIterator = boost::reverse_iterator<ConstIterator>;
    using Range = SimpleRandomAccessRangeFromPoint< ConstIterator, Iterator, DistanceFunctor >;
    using ConstRange = SimpleRandomAccessConstRangeFromPoint< ConstIterator, DistanceFunctor >;
    using Difference = std::ptrdiff_t;


    /// Policies as friend.
    friend TDomainPolicy<Self, Domain>;

    using TDomainPolicy<Self, Domain>::domain;
    
    /// Constructor.
    ImageView( MultiImage& aMultiImage, Label aLabel ) 
      : myMultiImage{&aMultiImage}, myLabel{aLabel}
    {}

    /// Default constructor.
    ImageView()
      : myMultiImage{nullptr}, myLabel{0}
    {}


    /// Get a value
    inline
    Value getValue( Point const& aPoint ) const
      {
        return myMultiImage->getValue( aPoint, myLabel );
      }

    /// Set a value
    inline
    void setValue( Point const& aPoint, Value aValue )
      {
        myMultiImage->setValue( aPoint, myLabel, aValue );
      }

    /// Get a value
    inline
    Value operator() ( Point const& aPoint ) const
      {
        return getValue( aPoint );
      }

    /**
     * Conversion to a ImageContainerBySTLVector
     * \todo is that the efficient way ? What about move syntax ?
     */
    explicit operator ImageContainerBySTLVector<Domain, Value>  () const
      {
        ImageContainerBySTLVector<Domain, Value> image{domain()};
        for ( auto const& point : domain() )
          image.setValue( point, getValue( point ) );
      
        return image;
      }

    inline
    Iterator begin()
      {
        return Iterator{this};
      }

    inline
    Iterator end()
      {
        return Iterator{this, true};
      }

    inline
    ConstIterator cbegin() const
      {
        return ConstIterator{this};
      }

    inline
    ConstIterator cend() const
      {
        return ConstIterator{this, true};
      }

    inline
    ReverseIterator rbegin()
      {
        return ReverseIterator{Iterator{this, true}};
      }

    inline
    ReverseIterator rend()
      {
        return ReverseIterator{Iterator{this}};
      }

    inline
    ConstReverseIterator crbegin() const
      {
        return ConstReverseIterator{ConstIterator{this, true}};
      }

    inline
    ConstReverseIterator crend() const
      {
        return ConstReverseIterator{ConstIterator{this}};
      }

    inline
    Range range()
      {
        return { begin(), end(), DistanceFunctor{cbegin()} };
      }

    inline
    ConstRange constRange() const
      {
        return { cbegin(), cend(), DistanceFunctor{cbegin()} };
      }

    /**
     * Iterator (not interoperable)
     */
    template < typename TImageView >
    class MyIterator
      : public boost::iterator_facade<
          MyIterator<TImageView>,
          Value,
          //boost::random_access_traversal_tag,
          std::random_access_iterator_tag,
          typename std::conditional< std::is_const<TImageView>::value, Value, typename MultiImage::Reference >::type
        >
      {
      public:
        /// Default constructor.
        MyIterator() noexcept {}

        /// Iterator from any point.
        MyIterator( TImageView* anImageView, Point const& aPoint ) noexcept
          : myImageView{ anImageView }
          {
            Domain domain = myImageView->domain();
            myExtent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
            myDimIndex = aPoint - domain.lowerBound();
            
            domain = myImageView->myMultiImage->domain();
            myGlobalExtent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
            myGlobalIndex = myImageView->myMultiImage->linearized( aPoint );

          }

        /// Copy constructor.
        template < typename TOtherImageView >
        MyIterator( MyIterator<TOtherImageView> const& other, 
            typename std::enable_if< std::is_convertible<TOtherImageView*, TImageView*>::value >::type* = 0 ) noexcept
          : myImageView{ other.myImageView }
          , myExtent{ other.myExtent }, myGlobalExtent{ other.myGlobalExtent }
          , myGlobalIndex{ other.myGlobalIndex }
          , myDimIndex{ other.myDimIndex }
          {}
        
        template < typename TOtherImageView >
        MyIterator( MyIterator<TOtherImageView> && other,
            typename std::enable_if< std::is_convertible<TOtherImageView*, TImageView*>::value >::type* = 0 ) noexcept
          : myImageView{ std::move(other.myImageView) }
          , myExtent{ std::move(other.myExtent) }, myGlobalExtent{ std::move(other.myGlobalExtent) }
          , myGlobalIndex{ std::move(other.myGlobalIndex) }
          , myDimIndex{ std::move(other.myDimIndex) }
          {}

        /// begin iterator.
        explicit MyIterator( TImageView* anImageView ) 
          : MyIterator{ anImageView, anImageView->domain().lowerBound() }
          {
          }

        /// end iterator.
        MyIterator( TImageView* anImageView, bool /* last */ ) 
          : MyIterator{ anImageView, anImageView->domain().upperBound() }
          {
            increment();
          }

      private:

        /// Interoperability
        template <class> friend class MyIterator;
        friend class boost::iterator_core_access;
        friend class DistanceFunctor;


        /// Increment of one step.
        void increment()
          {
            ++myGlobalIndex;
            ++myDimIndex[0];
            for ( std::size_t i = 1; i < Domain::dimension && myDimIndex[i-1] >= myExtent[i-1]; ++i )
              {
                myDimIndex[i-1] = 0;
                ++myDimIndex[i];
                std::size_t cum = myGlobalExtent[i-1] - myExtent[i-1];
                for ( std::size_t j = 0; j < i-1; ++j )
                  cum *= myGlobalExtent[j];

                myGlobalIndex += cum;
              }
          }

        /// Decrement of one step.
        void decrement()
          {
            --myGlobalIndex; /// \remarks Doesn't work if point coordinate are not signed ...
            --myDimIndex[0];
            for ( std::size_t i = 1; i < Domain::dimension && myDimIndex[i-1] < 0; ++i )
              {
                myDimIndex[i-1] = myExtent[i-1]-1;
                --myDimIndex[i];
                std::size_t cum = myGlobalExtent[i-1] - myExtent[i-1];
                for ( std::size_t j = 0; j < i-1; ++j )
                  cum *= myGlobalExtent[j];

                myGlobalIndex -= cum;
              }
          }

        template < typename TOtherImageView >
        inline
        bool equal( MyIterator<TOtherImageView> const& other ) const
          {
            return myGlobalIndex == other.myGlobalIndex;
          }

        template < class T = TImageView >
        inline
        typename std::enable_if< std::is_const<T>::value, Value >::type
        dereference() const
          {
            return myImageView->myMultiImage->getValueByIndex( myGlobalIndex, myImageView->myLabel );
          }

        template < class T = TImageView >
        inline
        typename std::enable_if< ! std::is_const<T>::value, Reference >::type
        dereference() const
          {
            return Reference( *(myImageView->myMultiImage), myImageView->domain().lowerBound()+myDimIndex, myImageView->myLabel, myGlobalIndex );
          }
        
        /// Distance to other iterator.
        std::ptrdiff_t distance_to( ConstIterator const& other ) const
          {
            const Point diff = other.myDimIndex - myDimIndex;
            std::ptrdiff_t dist = 0;
            for ( std::size_t i = Domain::dimension; i > 0 ; --i )
              dist = diff[i-1] + myExtent[i-1]*dist;

            return dist;
          }

        /// Distance to a point.
        std::ptrdiff_t distance_to( Point const& aPoint ) const
          {
            const Point diff = aPoint - myImageView->domain().lowerBound() - myDimIndex;
            std::ptrdiff_t dist = 0;
            for ( std::size_t i = Domain::dimension; i > 0 ; --i )
              dist = diff[i-1] + myExtent[i-1]*dist;

            return dist;
          }

        /// Advance by n steps. Not very efficient implementation ...
        void advance( std::ptrdiff_t n )
          {
            std::ptrdiff_t pos = 0;
            for ( std::size_t i = Domain::dimension; i > 0; --i )
              pos = myDimIndex[i-1] + myExtent[i-1]*pos;

            pos += n;

            for ( std::size_t i = 0; i < Domain::dimension; ++i )
              {
                myDimIndex[i] = pos % myExtent[i];
                pos /= myExtent[i];
              }

            myGlobalIndex = myImageView->myMultiImage->linearized( myImageView->domain().lowerBound() + myDimIndex );
          }

      private:
        TImageView* myImageView;
        Point myExtent, myGlobalExtent;
        typename Point::Coordinate myGlobalIndex;
        Point myDimIndex;
      };

    class DistanceFunctor
      {
      public:
        using Point = Self::Point;
        using Difference = Self::Difference;

        DistanceFunctor( ConstIterator const& anIterator ) noexcept
          : myIterator(anIterator)
          {}

        Difference operator() ( Point const& aPoint ) const noexcept
          {
            return myIterator.distance_to(aPoint);
          }

      private:
        const ConstIterator  myIterator;
      };

  private:
    MultiImage* myMultiImage;
    Label       myLabel;

  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

