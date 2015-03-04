#pragma once

#include <boost/concept_check.hpp>

namespace DGtal
{

namespace concepts
{

/*! Concept of value approximation
 */
template < typename TApprox >
struct CValueApproximation
  {
    typedef typename TApprox::value_type value_type;
    
    BOOST_CONCEPT_USAGE(CValueApproximation)
      {
        TApprox * approx = nullptr;
        value_type const* value = nullptr;
        bool test = approx->eval(*value);
        value_type const& default_value = approx->default_value;

        boost::ignore_unused_variable_warning(default_value);
        boost::ignore_unused_variable_warning(test);
      }
  };

} // namespace concepts


namespace approximation {

//! Base approximation class
template < typename T >
struct BaseValueApproximation
  {
    using value_type = T;
  };

//! No approximation with 0 default value
template < typename T >
struct NoValueApproximation : BaseValueApproximation<T>
  {
    static inline constexpr 
    bool eval ( T const& /* value */ ) { return false; }

    static const T default_value;
  };

template < typename T >
const T NoValueApproximation<T>::default_value = T(0);


//! Zero exact approximation
template < typename T >
struct ZeroValueApproximation : BaseValueApproximation<T>
  {
    using value_type = T;

    static inline constexpr 
    bool eval ( T const& value ) { return value == T(0); }

    static const T default_value;
  };

template < typename T > 
const T ZeroValueApproximation<T>::default_value = T(0);

//! Zero approximation given a tolerance
template < typename T >
struct ZeroTolValueApproximation : BaseValueApproximation<T>
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
struct NegativeValueApproximation : BaseValueApproximation<T>
  {
    static inline constexpr 
    bool eval ( T const& value ) { return value <= T(0); }
    
    static const T default_value;
  };

template < typename T >
const T NegativeValueApproximation<T>::default_value = T(0);

//! Negative value approximation given a tolerance ( value <= ? tol )
template < typename T >
struct NegativeTolValueApproximation : BaseValueApproximation<T>
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

