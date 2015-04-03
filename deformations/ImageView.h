#pragma once

#include <cstddef>
#include <array>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

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
        return static_cast<TImageView const*>(this)->myMultiImage.domain();
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
        return image.myMultiImage.getBoundingBox( image.myLabel, myBuffer );
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

    // Iterators & Ranges
    class ConstIterator;
    using Iterator        = ConstIterator;
    using ReverseIterator = boost::reverse_iterator<Iterator>;
    using ConstReverseIterator = boost::reverse_iterator<ConstIterator>;
    class ConstRange;
    using Range  = ConstRange;


    /// Policies as friend.
    friend TDomainPolicy<Self, Domain>;

    using TDomainPolicy<Self, Domain>::domain;
    
    /// Constructor.
    ImageView( MultiImage& aMultiImage, Label aLabel ) 
      : myMultiImage{aMultiImage}, myLabel{aLabel}
    {}

    /// No default constructor.
    ImageView() = delete;

    /// Get a value
    inline
    Value getValue( Point const& aPoint ) const
      {
        return myMultiImage.getValue( aPoint, myLabel );
      }

    /// Set a value
    inline
    void setValue( Point const& aPoint, Value aValue )
      {
        myMultiImage.setValue( aPoint, myLabel, aValue );
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
    Iterator begin() const
      {
        return ConstIterator{this};
      }

    inline
    Iterator end() const
      {
        return ++ConstIterator{this, true};
      }

    inline
    ConstIterator cbegin() const
      {
        return ConstIterator{this};
      }

    inline
    ConstIterator cend() const
      {
        return ++ConstIterator{this, true};
      }

    inline
    ReverseIterator rbegin() const
      {
        return ReverseIterator{ConstIterator{this, true}};
      }

    inline
    ConstReverseIterator rend() const
      {
        return ReverseIterator{--ConstIterator{this}};
      }

    inline
    Range range() const
      {
        return { begin(), end() };
      }

    inline
    ConstRange constRange() const
      {
        return { cbegin(), cend() };
      }

    /**
     * Constant Iterator
     */
    class ConstIterator
      : public boost::iterator_facade<
          ConstIterator,
          Value,
          boost::bidirectional_traversal_tag,
          Value
        >
      {
      public:
        ConstIterator() {}

        explicit ConstIterator( Self const* anImageView ) 
          : myImageView{ anImageView }
          {
            Domain domain = myImageView->domain();
            myExtent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
            myGlobalIndex = myImageView->myMultiImage.linearized( domain.lowerBound() );
            
            domain = myImageView->myMultiImage.domain();
            myGlobalExtent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);

            myDimIndex = Point::diagonal(0);
          }

        ConstIterator( Self const* anImageView, bool /* last */ ) 
          : myImageView{ anImageView }
          {
            Domain domain = myImageView->domain();
            myExtent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
            myGlobalIndex = myImageView->myMultiImage.linearized( domain.upperBound() );
            
            domain = myImageView->myMultiImage.domain();
            myGlobalExtent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);

            myDimIndex = myExtent - Point::diagonal(1);
          }

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

        void decrement()
          {
            --myGlobalIndex; /// \remarks Doesn't is point coordinate are not signed ...
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

        inline
        bool equal( ConstIterator const& other ) const
          {
            return myGlobalIndex == other.myGlobalIndex;
          }

        inline
        Value dereference() const
          {
            return myImageView->myMultiImage.getValueByIndex( myGlobalIndex, myImageView->myLabel );
          }

      private:
        Self const* myImageView;
        Point myExtent, myGlobalExtent;
        typename Point::Coordinate myGlobalIndex;
        Point myDimIndex;
      };

    class ConstRange
      : public boost::iterator_range<ConstIterator>
      {
      public:
        using ConstIterator = ConstIterator;
      };

  private:
    MultiImage& myMultiImage;
    Label       myLabel;

  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

