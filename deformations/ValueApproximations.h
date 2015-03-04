#pragma once

namespace DGtal
{

namespace approximation {

//! No approximation with 0 default value
template < typename T >
struct NoValueApproximation
  {
    static inline constexpr 
    bool eval ( T const& /* value */ ) { return false; }

    static const T default_value;
  };

template < typename T >
const T NoValueApproximation<T>::default_value = T(0);


//! Zero exact approximation
template < typename T >
struct ZeroValueApproximation
  {
    static inline constexpr 
    bool eval ( T const& value ) { return value == T(0); }

    static const T default_value;
  };

template < typename T > 
const T ZeroValueApproximation<T>::default_value = T(0);

//! Zero approximation given a tolerance
template < typename T >
struct ZeroTolValueApproximation
  {
    constexpr ZeroTolValueApproximation( T const& zero_tol ) : tol{zero_tol} {}
    inline constexpr 
    bool eval (T const& value ) const { return -tol <= value && value <= tol; }

    T tol;
    static const T default_value;
  };

template < typename T >
const T ZeroTolValueApproximation<T>::default_value = T(0);

//! Negative value approximation
template < typename T >
struct NegativeValueApproximation
  {
    static inline constexpr 
    bool eval ( T const& value ) { return value <= T(0); }
    
    static const T default_value;
  };

template < typename T >
const T NegativeValueApproximation<T>::default_value = T(0);

//! Negative value approximation given a tolerance ( value <= ? tol )
template < typename T >
struct NegativeTolValueApproximation
  {
    constexpr NegativeTolValueApproximation( T const& zero_tol ) : tol{zero_tol} {}
    inline constexpr 
    bool eval ( T const& value ) const { return value <= tol; }

    T tol;
    static const T default_value;
  };

template < typename T >
const T NegativeTolValueApproximation<T>::default_value = T(0);

} // namespace approximation

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

