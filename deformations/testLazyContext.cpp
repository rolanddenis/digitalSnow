#include <iostream>
#include <string>
#include <utility>
#include <tuple>
#include <cmath>

template < typename T, std::size_t I >
struct TypeAlias
{
  const T value;
  operator T() const { return value; }
  T const& operator() () const { return value; }
};

struct ContextTraits
{
  using A = TypeAlias<double, 0>;
  using B = TypeAlias<int, 1>;
  using C = TypeAlias<double, 2>;
  using D = TypeAlias<std::string, 3>;
};

struct Context
{
  double a;
  int b;
  double c;
  std::string d;
};

/*
template <
  typename AFunctor,
  typename BFunctor,
  typename CFunctor,
  typename DFunctor
>
struct LazyContext
{
  AFunctor a;
  BFunctor b;
  CFunctor c;
  DFunctor d;

  operator ContextTraits::A() const { return a(); }
  operator ContextTraits::B() const { return b(); }
  operator ContextTraits::C() const { return c(); }
  operator ContextTraits::D() const { return d(); }
  operator Context() const { return { a(), b(), c(), d() }; }
};


template <
  typename AFunctor,
  typename BFunctor,
  typename CFunctor,
  typename DFunctor
>
LazyContext<AFunctor, BFunctor, CFunctor, DFunctor>
makeLazyContext( AFunctor && a, BFunctor && b, CFunctor && c, DFunctor && d )
{
  return {
      std::forward<AFunctor>(a),
      std::forward<BFunctor>(b),
      std::forward<CFunctor>(c),
      std::forward<DFunctor>(d)
  };
}
*/

template < typename T > struct FunctorTraitsImpl;

template < typename Class, typename Ret, typename... Args >
struct FunctorTraitsImpl< Ret (Class::*) (Args...) const >
{
  using Arguments = std::tuple< typename std::decay<Args>::type... >;
};

template < typename T, bool IsAFunction >
struct FunctionTraitsImpl
  : FunctorTraitsImpl< decltype(&T::operator()) >
{};

template < typename Ret, typename... Args >
struct FunctionTraitsImpl< Ret(Args...), true >
{
  using Arguments = std::tuple< typename std::decay<Args>::type... >;
};


template < typename T >
struct FunctionTraits
  : FunctionTraitsImpl< T, std::is_function<T>::value >
{};

template < typename Functor >
struct Operator
{
  Functor aFunctor;

  // Context calculators
  ContextTraits::A  calcContext( double seed, ContextTraits::A const* ) const { return { std::cos(seed) }; }
  ContextTraits::B  calcContext( double seed, ContextTraits::B const* ) const { return { static_cast<int>(10*std::sin(seed)) }; }
  ContextTraits::C  calcContext( double seed, ContextTraits::C const* ) const { return { std::tan(seed) }; }
  ContextTraits::D  calcContext( double seed, ContextTraits::D const* ) const { return { std::to_string(seed) }; }

  template < typename... Args >
  double applyFunctorWithContext( double seed, std::tuple<Args...> const* ) const
    {
      return aFunctor( calcContext( seed, static_cast<Args const*>(nullptr) )... );
    }

  double getValue( double seed ) const
    {
      return applyFunctorWithContext( seed, static_cast<typename FunctionTraits<Functor>::Arguments const*>(nullptr) );
    }
};

template < typename Functor >
Operator<Functor>
makeOperator( Functor && aFunctor )
{
  return { std::forward<Functor>(aFunctor) };
}


int main( int argc, char* argv[] )
{
  const double seed = std::stod( argv[1] );
  
  auto const fn = makeOperator( [] ( ContextTraits::B const& , ContextTraits::A const& a ) mutable { return a(); } );

  std::cout << fn.getValue( seed ) << std::endl;

  return 0;
}
