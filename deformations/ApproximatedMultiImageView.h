#pragma once

#include <cstddef>
#include <array>
#include <type_traits>

#include <DGtal/images/ImageContainerBySTLVector.h> // For conversion purpose
#include <DGtal/images/ArrayImageIterator.h>
#include <DGtal/base/IteratorCompletion.h>
#include <DGtal/kernel/domains/Linearizer.h>

namespace DGtal
{

namespace image_view
{

/**
 * Policy that returns the full domain
 */
template < typename TApproximatedMultiImageView, typename TDomain >
struct FullDomain
  {
    inline
    TDomain domain() const
      {
        return static_cast<TApproximatedMultiImageView const*>(this)->myMultiImage->domain();
      }
  protected:
    ~FullDomain() {}
  };

/**
 * Policy that returns the bounding box (with buffer) as the image domain.
 */
template < typename TApproximatedMultiImageView, typename TDomain >
class BoundingBoxAsDomain
  {
  public:
    using Point = typename TDomain::Point;

    inline
    TDomain domain() const
      {
        TApproximatedMultiImageView const& image = * static_cast<TApproximatedMultiImageView const*>(this);
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
  
  protected:
    ~BoundingBoxAsDomain() {}
  };

} // namespace image_view


/**
 * Image View for ApproximatedMultiImage (but not only ?)
 */
template <
  typename TMultiImage,
  template<typename,typename> class TDomainPolicy = image_view::FullDomain
>
class ApproximatedMultiImageView
    : public TDomainPolicy< ApproximatedMultiImageView<TMultiImage, TDomainPolicy>, typename TMultiImage::Domain >
    , public IteratorCompletion< ApproximatedMultiImageView<TMultiImage, TDomainPolicy> >
  {
  public:
    // Typedefs
    using Self        = ApproximatedMultiImageView<TMultiImage, TDomainPolicy>;
    using MultiImage  = TMultiImage;
    using Point       = typename MultiImage::Point;
    using Label       = typename MultiImage::Label;
    using Domain      = typename MultiImage::Domain;
    using Value       = typename MultiImage::Value;
    using Reference   = typename MultiImage::Reference;

    // Iterators & Ranges
    template<class> friend class ArrayImageIterator;
    using Iterator = typename IteratorCompletionTraits<Self>::Iterator;
    using ConstIterator = typename IteratorCompletionTraits<Self>::ConstIterator;
    
    /// Policies as friend.
    friend TDomainPolicy<Self, Domain>;
    using TDomainPolicy<Self, Domain>::domain;
    
    /// Constructor.
    ApproximatedMultiImageView( MultiImage& aMultiImage, Label aLabel ) 
      : myMultiImage{&aMultiImage}, myLabel{aLabel}
    {}

    /// Default constructor.
   ApproximatedMultiImageView()
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
    Value setValue( Point const& aPoint, Value aValue )
      {
        return myMultiImage->setValue( aPoint, myLabel, aValue );
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
        return Iterator{this, myMultiImage->domain(), domain()};
      }

    inline
    Iterator end()
      {
        return Iterator{this, myMultiImage->domain(), domain(), true};
      }
    
    inline
    ConstIterator begin() const
      {
        return ConstIterator{this, myMultiImage->domain(), domain()};
      }

    inline
    ConstIterator end() const
      {
        return ConstIterator{this, myMultiImage->domain(), domain(), true};
      }

    inline
    ConstIterator cbegin() const
      {
        return ConstIterator{this, myMultiImage->domain(), domain()};
      }

    inline
    ConstIterator cend() const
      {
        return ConstIterator{this, myMultiImage->domain(), domain(), true};
      }

  public: // Should be private since ArrayImageIterator is a friend but g++ 4.9.1 don't care ... (no prob with clang++)
    
    inline
    Value dereference( Point const& /* aPoint */ , typename Point::Coordinate aFullIndex ) const
      {
        return myMultiImage->getValueByIndex( aFullIndex, myLabel );
      }

    inline
    Reference dereference( Point const& aPoint, typename Point::Coordinate aFullIndex )
      {
        return Reference( *myMultiImage, aPoint, myLabel, aFullIndex );
      }
    
  private:
    MultiImage* myMultiImage;
    Label       myLabel;

  };

template <
  typename TMultiImage,
  template<typename,typename> class TDomainPolicy
>
class IteratorCompletionTraits< ApproximatedMultiImageView<TMultiImage, TDomainPolicy> >
  {
public:
    using Self          = ApproximatedMultiImageView<TMultiImage, TDomainPolicy>;
    using ConstIterator = ArrayImageIterator<Self const>;
    using Iterator      = typename std::conditional< std::is_const<TMultiImage>::value, ConstIterator, ArrayImageIterator<Self> >::type;
    
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

