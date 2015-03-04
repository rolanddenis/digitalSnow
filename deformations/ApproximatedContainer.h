#pragma once

#include <boost/concept_check.hpp>
#include <cstddef>
#include <type_traits>

#include "ValueApproximations.h"

namespace DGtal
{

namespace concepts
{
  // TODO
} // namespace concepts

template <
    typename TContainer,
    typename TApproximation = approximations::NoValueApproximation< typename TContainer::value_type >
>
class ApproximatedContainer
  {

    BOOST_CONCEPT_ASSERT( (concepts::CValueApproximation<TApproximation>) );
    static_assert( std::is_same< typename TApproximation::value_type, typename TContainer::value_type >::value, "Container and approximation must have same value_type");

  public:
    using container_type      = TContainer;
    using approximation_type  = TApproximation;
    using size_type           = typename container_type::size_type;
    using value_type          = typename container_type::value_type;
    using reference           = typename container_type::reference;
    using const_reference     = typename container_type::const_reference;
    
    ApproximatedContainer(size_type n = 0, approximation_type const& approx = approximation_type() ) : m_approx{approx}, m_data(0)
      {
        resize(n);
      }

    inline
    void resize( size_type n, value_type val = m_approx.default_value )
      {
        m_size = n;
        m_data.resize(m_size, val);
      }

    inline 
    size_type size() const
      {
        return m_size;
      }

    inline
    reference operator[] (size_type i)
      {
        return m_data[i];
      }

    inline
    const_reference operator[] (size_type i) const
      {
        return m_data[i];
      }

  private:
    approximation_type m_approx;
    container_type m_data;
    size_type m_size;

  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

