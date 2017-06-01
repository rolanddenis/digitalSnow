template <typename TDerived>
class ImageOperatorTraits;

template <typename TDerived>
class ImageOperator
{
  using Derived = TDerived;

  ~ImageOperator = delete;

  Derived      & getDerived()       { return *static_cast<Derived*>(this);       }
  Derived const& getDerived() const { return *static_cast<Derived const*>(this); }

  template < typename TImage >
  ... operator() ( TImage const& anImage ) { 
};
