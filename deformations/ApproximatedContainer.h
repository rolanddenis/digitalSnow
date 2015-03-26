#pragma once

#include <boost/concept_check.hpp>
#include <cstddef>
#include <type_traits>

#include "ValueApproximations.h"
#include <DGtal/base/LabelledMap.h>

#include <iostream>

namespace DGtal
{

namespace concepts
{
  // TODO
} // namespace concepts


template <
    typename TContainer,
    typename TApproximation = approximations::NoValueApproximation< typename TContainer::value_type > ///@warning Incompatible with LabelledMap & Map !!!
>
class ApproximatedContainer;
/*
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
    void resize( size_type n, value_type val )
      {
        m_size = n;
        m_data.resize(m_size, val);
      }

    inline
    void resize( size_type n )
      {
        resize( n, m_approx.default_value );
      }

    inline 
    size_type size() const
      {
        return m_size;
      }

    inline
    size_type stored_size() const
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
*/



template <
    typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M,
    typename TApproximation
>
class ApproximatedContainer< DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation >
  {

    BOOST_CONCEPT_ASSERT( (concepts::CValueApproximation<TApproximation>) );
    static_assert( std::is_same< typename TApproximation::value_type, TData>::value, "Container and approximation must have same value_type");

  public:
    template <typename> class ApproximatedReference;

    using self_type           = ApproximatedContainer< DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation >;
    using container_type      = DGtal::LabelledMap<TData, L, TWord, N, M>;
    using approximation_type  = TApproximation;
    using size_type           = typename container_type::size_type;
    using value_type          = typename container_type::mapped_type;
    using reference           = ApproximatedReference<self_type>;
    using const_reference     = value_type const&;
    
    ApproximatedContainer(size_type n = 0, approximation_type const& approx = approximation_type() ) : m_approx{approx}, m_data{}, m_size(0)
      {
        resize(n);
      }

    void resize( size_type n, value_type val )
      {
        if ( n < m_size )
          {
            m_data.erase( m_data.lower_bound(n), m_data.end() );
          }
        else if (n > m_size && val != m_approx.default_value)
          {
              for ( size_type i = m_size ; i < n ; ++i )  m_data[i] = val;
          }

        m_size = n;
      }

    inline
    void resize( size_type n )
      {
        resize( n, m_approx.default_value );
      }

    inline 
    size_type size() const noexcept
      {
        return m_size;
      }

    inline
    size_type stored_size() const noexcept
      {
        return m_data.size();
      }

    inline
    reference operator[] (size_type i)
      {
        return reference{*this, i};
      }

    inline
    const_reference operator[] (size_type i) const noexcept
      {
        //return reference{*static_cast<const self_type*>(this),i};
        //return ApproximatedReference<const self_type>{*this, i};
        if ( m_data.count(i) > 0 )
          {
            return m_data.fastAt(i);
          }
        else
          {
            return m_approx.default_value;
          }
      }

    /**
     *
     * \remark Missing mutable iterator for LabelledMap.
     */
    /*
    class ApproximatedLazyReference
      {
      public:
        ApproximatedLazyReference( size_type i )
          {
            const auto it = m_data.find(i);
            if ( it != m_data.end() )
              {
                m_exist = true;
                m_val

      private:
        bool m_exist;
        value_type m_value;
        typename container_type::Value find

      };
    */

    template <typename TApproximatedContainer>
    class ApproximatedReference
      {
      public:
        using self_type = ApproximatedReference<TApproximatedContainer>; ///< Type of this class.
        using value_type = typename TApproximatedContainer::value_type;  ///< Type of stored values.

        /// Constructor
        ApproximatedReference( TApproximatedContainer & container, size_type i ) : m_container{container}, m_index{i} {}

        /// Implicit conversion to the underlying value.
        operator value_type () const
          {
            const auto it = m_container.m_data.find(m_index);
            if ( it != m_container.m_data.end() )
              {
                return (*it).second;
              }
            else
              {
                return m_container.m_approx.default_value;
              }
          }

        /// Modification of the associated value.
        self_type& operator= ( value_type value )
          {
            if ( m_container.m_approx.eval(value) )
              {
                const auto it = m_container.m_data.find(m_index);
                if ( it != m_container.m_data.end() )
                  {
                    m_container.m_data.erase(it);
                  }
              }
            else
              {
                m_container.m_data[m_index] = value;
              }
            
            return *this;
          }


      private:
        TApproximatedContainer & m_container; ///< Associated container
        const size_type m_index; ///< Index of the referenced value in the container.
      };

  private:
    approximation_type m_approx;
    container_type m_data;
    size_type m_size;

  };


















} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

