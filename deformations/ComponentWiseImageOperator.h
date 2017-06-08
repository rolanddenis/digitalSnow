#include "ImageOperator.h"

template < typename TDerived, typename TImage >
struct ComponentWiseImageContext
{
  using Derived = TDerived;
  using Image   = TImage;
  using Value   = typename Image::Value;
  using Point   = typename Image::Point;
  
  Derived const& getDerived() const& { return *static_cast<Derived const*>(this); }

  const Point point() const { return getDerived().point(); }
  const Value value() const { return getDerived().value(); }
};

template < typename TOperatorStorage, typename TImage >
class ComponentWiseImageOperatorResult;

template < typename TOperatorStorage, typename TImage >
struct ImageOperatorResultTraits< ComponentWiseImageOperatorResult<TOperatorStorage, TImage > >
{
  using Image = TImage;
};

template < typename TOperatorStorage, typename TImage >
class ComponentWiseImageOperatorResult
  : public ImageOperatorResult< ComponentWiseImageOperatorResult<TOperatorStorage, TImage> >
{
public:
  using Image    = TImage;
  using Point    = typename Image::Point;
  using OperatorStorage = TOperatorStorage;
  using Operator = typename std::decay<OperatorStorage>::type;
  using Context  = ComponentWiseImageContext<Point, double>;

private:
  OperatorStorage myOperator;
  Image    const& myImage;

public:
  ComponentWiseImageOperatorResult( Operator const& anOperator, Image const& anImage )
    : myOperator( anOperator ), myImage( anImage )
    {
    }

  double getValue( typename Image::Point const& aPoint ) const
    {
      //return myOperator.applyOnImageContext( Context{ aPoint, myImage(aPoint) } );
      return myOperator.applyOnImageContext( Context{ aPoint + Point::diagonal(static_cast<int>(std::tan(aPoint[0]))), myImage(aPoint) } );
      //return myOperator.myFunctor( myImage(aPoint) );

    }

};

template < typename TFunctorStorage >
class ComponentWiseImageOperator;

template < typename TFunctorStorage >
struct ImageOperatorTraits< ComponentWiseImageOperator<TFunctorStorage> >
{
  using Operator     = ComponentWiseImageOperator<TFunctorStorage>;

  template < typename TImage >
  using RValueResult = ComponentWiseImageOperatorResult< Operator const , typename std::decay<TImage>::type >;
  
  template < typename TImage >
  using LValueResult = ComponentWiseImageOperatorResult< Operator const&, typename std::decay<TImage>::type >;
};



template < typename TFunctorStorage >
class ComponentWiseImageOperator
  : public ImageOperator< ComponentWiseImageOperator<TFunctorStorage> >
{
public:
  using FunctorStorage = TFunctorStorage;
  using Functor = typename std::decay<FunctorStorage>::type;
  using Self    = ComponentWiseImageOperator<FunctorStorage>;

private:
  FunctorStorage myFunctor;

public:
  template < typename TFunctor2 >
  ComponentWiseImageOperator( TFunctor2 && aFunctor, int )
    : myFunctor( std::forward<TFunctor2>(aFunctor) )
    {}

  template < typename TImage >
  typename ImageOperatorTraits<Self>::template LValueResult<TImage>
  applyOnImage ( TImage && anImage ) const&
    {
      return { *this, std::forward<TImage>( anImage ) };
    }
  
  template < typename TImage >
  typename ImageOperatorTraits<Self>::template RValueResult<TImage>
  applyOnImage ( TImage && anImage ) &&
    {
      return { std::move(*this), std::forward<TImage>( anImage ) };
    }

public:
  template < typename TImageContext >
  auto applyOnImageContext( TImageContext const& aContext ) const
      -> decltype( myFunctor( aContext ) )
    {
      return myFunctor( aContext );
    }
};


template < typename TFunctor >
ComponentWiseImageOperator< TFunctor >
makeComponentWiseImageOperator ( TFunctor && aFunctor )
{
  return { std::forward<TFunctor>(aFunctor), 0 };
}
