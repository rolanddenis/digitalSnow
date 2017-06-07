template < typename TImageOperator >
struct ImageOperatorTraits;

template <typename TDerived>
class ImageOperator
{
public:
  using Derived = TDerived;

protected:
  ~ImageOperator() = default;

public:
  Derived     && getDerived()     && { return static_cast<Derived&&>(*this);       }
  Derived const& getDerived() const& { return *static_cast<Derived const*>(this); }

  template < typename TImage >
  typename ImageOperatorTraits<Derived>::template LValueResult<TImage>
  operator() ( TImage && anImage ) const&
    {
      return getDerived().applyOnImage( std::forward<TImage>(anImage) );
    }
  
  template < typename TImage >
  typename ImageOperatorTraits<Derived>::template RValueResult<TImage>
  operator() ( TImage && anImage ) &&
    {
      return std::move(*this).getDerived().applyOnImage( std::forward<TImage>(anImage) );
    }
};

template < typename TImageOperatorResult >
struct ImageOperatorResultTraits;

template <typename TDerived>
class ImageOperatorResult
{
public:
  using Derived = TDerived;
  using Image   = typename ImageOperatorResultTraits<Derived>::Image;

protected:
  ~ImageOperatorResult() = default;

public:
  Derived const& getDerived() const { return *static_cast<Derived const*>(this); }

  double operator() ( typename Image::Point const& aPoint ) const
    {
      return getDerived().getValue( aPoint );
    }
};
