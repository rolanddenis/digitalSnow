#include "ImageOperator.h"

template < typename TPoint, typename TValue >
struct ComponentWiseImageContext
{
  using Value = TValue;
  using Point = TPoint;

  Point point;
  Value value;
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
      return myOperator.applyOnImageContext( Context{ aPoint + Point::diagonal(static_cast<int>(std::sin(aPoint[0]))), myImage(aPoint) } );
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

  //template < typename TOperatorStorage, typename TImage >
  //friend class ComponentWiseImageOperatorResult;

private:
  FunctorStorage myFunctor;

public:
  ComponentWiseImageOperator( Functor const& aFunctor )
    : myFunctor( aFunctor )
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
ComponentWiseImageOperator< TFunctor const& >
makeComponentWiseImageOperator ( TFunctor const& aFunctor )
{
  return { aFunctor };
}

template < typename TFunctor >
ComponentWiseImageOperator< TFunctor const >
makeComponentWiseImageOperator ( TFunctor && aFunctor )
{
  return { std::forward<TFunctor>(aFunctor) };
}
