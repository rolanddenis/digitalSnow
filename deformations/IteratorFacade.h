#pragma once

#include <boost/iterator/reverse_iterator.hpp>
#include <DGtal/base/SimpleRandomAccessRangeFromPoint.h>
#include <DGtal/base/SimpleRandomAccessConstRangeFromPoint.h>

namespace DGtal
{

template < typename TDerived >
class IteratorTraits;

template < 
  typename TDerived
>
class IteratorFacade
  {
  public:

    using Iterator        = typename IteratorTraits<TDerived>::Iterator;
    using ConstIterator   = typename IteratorTraits<TDerived>::ConstIterator;
    using DistanceFunctor = typename IteratorTraits<TDerived>::DistanceFunctor;

    using ReverseIterator       = boost::reverse_iterator<Iterator>;
    using ConstReverseIterator  = boost::reverse_iterator<ConstIterator>;
    using Range = SimpleRandomAccessRangeFromPoint< ConstIterator, Iterator, DistanceFunctor >;
    using ConstRange = SimpleRandomAccessConstRangeFromPoint< ConstIterator, DistanceFunctor >;
    using Difference = std::ptrdiff_t;

    ReverseIterator rbegin()
      {
        return ReverseIterator{ static_cast<TDerived*>(this)->end() };
      }

    inline
    ReverseIterator rend()
      {
        return ReverseIterator{ static_cast<TDerived*>(this)->begin() };
      }

    inline
    ConstReverseIterator crbegin() const
      {
        return ConstReverseIterator{ static_cast<TDerived*>(this)->cend() };
      }

    inline
    ConstReverseIterator crend() const
      {
        return ConstReverseIterator{ static_cast<TDerived*>(this)->cbegin() };
      }

    inline
    Range range()
      {
        TDerived* const derived = static_cast<TDerived*>(this);
        return { 
            derived->begin(), 
            derived->end(),
            typename IteratorTraits<TDerived>::DistanceFunctor( derived )
        };
      }

    inline
    ConstRange constRange() const
      {
        TDerived const* const derived = static_cast<TDerived const*>(this);
        return { 
            derived->cbegin(), 
            derived->cend(), 
            typename IteratorTraits<TDerived>::DistanceFunctor( derived )
        };
      }

  protected:
    ~IteratorFacade()
      {}
  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

