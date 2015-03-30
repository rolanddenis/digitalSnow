#pragma once

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

  private:
    MultiImage& myMultiImage;
    Label       myLabel;

  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

