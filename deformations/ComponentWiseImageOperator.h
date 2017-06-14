#include "ImageOperator.h"
#include "FunctionTraits.h"
#include "ImageTransform.h"

/** Dummy type in order to identify arguments.
 *
 * @tparam  T   Underlying type.
 * @tparam  I   Argument unique identifier.
 */
template < typename T, std::size_t I >
struct TypeProxy
{
  const T value; ///< Underlying value.
  operator T() const { return value; }           ///< Underlying value accessor by evaluation.
  T const& operator() () const { return value; } ///< Underlying value accessor by implicit conversion.
};

/** Operator context for ComponentWiseImageOperator.
 *
 * @tparam TImage Image type.
 */
template < typename TImage >
struct ComponentWiseImageContext
{
  using Image   = TImage;
  using Value   = TypeProxy<typename Image::Value, 0>; ///< Unique type identifying the first argument (image's value).
  using Point   = TypeProxy<typename Image::Point, 1>; ///< Unique type identifying the second argument (space point).
  //using Id      = TypeProxy<std::size_t, 2>;           ///< Point id in the domain (position in a continuous storage).
};

template < typename TOperatorStorage, typename TImage >
class ComponentWiseImageOperatorResult;

template < typename TOperatorStorage, typename TImage >
struct ImageOperatorResultTraits< ComponentWiseImageOperatorResult<TOperatorStorage, TImage > >
{
  using Image = TImage;
};

/** Result proxy of ComponentWiseImageOperator.
 *
 * @tparam  TOperatorStorage    Storage type of the ComponentWiseImageOperator.
 * @tparam  TImage              Type of the image on which is applied the operator.
 */
template < typename TOperatorStorage, typename TImage >
class ComponentWiseImageOperatorResult
  : public ImageOperatorResult< ComponentWiseImageOperatorResult<TOperatorStorage, TImage> >
{
public:
  using Image    = TImage;
  using Value    = typename Image::Value;
  using Point    = typename Image::Point;
  using OperatorStorage = TOperatorStorage;
  using Operator = typename std::decay<OperatorStorage>::type;
  using Context  = ComponentWiseImageContext<Image>;

private:
  OperatorStorage myOperator;
  Image    const& myInputImage;

private:
  // Context calculators
  typename Context::Value calcContext( Point const& , Value const& aValue, typename Context::Value const* ) const { return { aValue }; } ///< Image value calculator.
  typename Context::Point calcContext( Point const& aPoint, Value const&,  typename Context::Point const* ) const { return { aPoint }; } ///< Space point calculator.

  template < typename... Args >
  double applyOperatorWithContext(
        Point const& aPoint, Value const& aValue,
        std::tuple<Args...> const* ) const
    {
      return myOperator.applyWithContext( calcContext( aPoint, aValue, static_cast<Args const*>(nullptr) )... );
    }

public:
  /** Constructor.
   * @param anOperator  Applied operator.
   * @param anImage     Image on which the operator is applied.
   */
  ComponentWiseImageOperatorResult( Operator const& anOperator, Image const& anInputImage )
    : myOperator( anOperator ), myInputImage( anInputImage )
    {
    }

  /// Evaluates the operator on the given point of the image.
  double getValue( Point const& aPoint ) const
    {
      return applyOperatorWithContext(
        aPoint, myInputImage(aPoint),
        static_cast<typename Operator::FunctorArgs const*>(nullptr) );
    }

  /// Saving result in an image
  /// Rmk: it should be in ImageOperatorResult base class since the behavior is
  ///      the same for almost all operators.
  /// TODO
  template < typename TOutputImage >
  void saveTo( TOutputImage & anOutputImage ) const
  {
    ASSERT( anOutputImage.domain() == myInputImage.domain() );

    ImageTransform( myInputImage, anOutputImage,
        [this] ( Point const& aPoint, Value const& aValue )
        {
            return applyOperatorWithContext(
                aPoint, aValue,
                static_cast<typename Operator::FunctorArgs const*>(nullptr) );
        }
    );
  }
};


template < typename TFunctorStorage >
class ComponentWiseImageOperatorClass;

/// Image operator traits for ComponentWiseImageOperator.
template < typename TFunctorStorage >
struct ImageOperatorTraits< ComponentWiseImageOperatorClass<TFunctorStorage> >
{
  using Operator     = ComponentWiseImageOperatorClass<TFunctorStorage>;

  /// Result type when the operator is a rvalue.
  template < typename TImage >
  using RValueResult = ComponentWiseImageOperatorResult< Operator const , typename std::decay<TImage>::type >;

  /// Result type when the operator is a lvalue.
  template < typename TImage >
  using LValueResult = ComponentWiseImageOperatorResult< Operator const&, typename std::decay<TImage>::type >;
};


/** Component-wise image operator.
 *
 * @tparam  TFunctorStorage Storage type of the underlying functor.
 */
template < typename TFunctorStorage >
class ComponentWiseImageOperatorClass
  : public ImageOperator< ComponentWiseImageOperatorClass<TFunctorStorage> >
{
public:
  using FunctorStorage = TFunctorStorage;
  using Functor     = typename std::decay<FunctorStorage>::type;
  using FunctorArgs = typename FunctionTraits<Functor>::Arguments; ///< Functors arguments types.
  using Self    = ComponentWiseImageOperatorClass<FunctorStorage>;

private:
  FunctorStorage myFunctor;

public:
  template < typename TFunctor2 >
  ComponentWiseImageOperatorClass( TFunctor2 && aFunctor, int )
    : myFunctor( std::forward<TFunctor2>(aFunctor) )
    {}

  /// Operator application on lvalue context.
  template < typename TImage >
  typename ImageOperatorTraits<Self>::template LValueResult<TImage>
  applyOnImage ( TImage && anImage ) const&
    {
      return { *this, std::forward<TImage>( anImage ) };
    }

  /// Operator application on rvalue context.
  template < typename TImage >
  typename ImageOperatorTraits<Self>::template RValueResult<TImage>
  applyOnImage ( TImage && anImage ) &&
    {
      return { std::move(*this), std::forward<TImage>( anImage ) };
    }

public:
  /// Apply underlying functor with given args.
  template < typename... TArgs >
  auto applyWithContext( TArgs const&... someArgs ) const
      -> decltype( myFunctor( someArgs... ) )
    {
      return myFunctor( someArgs... );
    }
};


/// ComponentWiseImageOperator construction helper.
template < typename TFunctor >
ComponentWiseImageOperatorClass< TFunctor >
ComponentWiseImageOperator ( TFunctor && aFunctor )
{
  return { std::forward<TFunctor>(aFunctor), 0 };
}

/// Image input in a component-wise image operator.
/// TODO: check image.
template < typename TImage, typename TFunctor >
auto operator >> ( TImage && anInputImage, ComponentWiseImageOperatorClass<TFunctor> const& anImageOperator )
    -> decltype( anImageOperator( std::forward<TImage>(anInputImage) ) )
{
    return anImageOperator( std::forward<TImage>(anInputImage) );
}

/// Image output from a component-wise image operator result.
/// TODO: check image.
template < typename TOperator, typename TInputImage, typename TOutputImage >
TOutputImage operator >> ( ComponentWiseImageOperatorResult<TOperator, TInputImage> const& anImageOperatorResult, TOutputImage & anOutputImage )
{
    anImageOperatorResult.saveTo( anOutputImage );
    return anOutputImage;
}
