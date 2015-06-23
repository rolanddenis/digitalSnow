#pragma once

#include <boost/iterator/reverse_iterator.hpp>
#include <DGtal/base/SimpleRandomAccessRangeFromPoint.h>
#include <DGtal/base/SimpleRandomAccessConstRangeFromPoint.h>

namespace DGtal
{

/** Traits that must be specialized for each IteratorFace derived class.
 *
 * This traits must show:
 * - a typedef Iterator corresponding to the derived class mutable iterator.
 * - a typedef ConstIterator corresponding to the derived class constant iterator.
 * - a class DistanceFunctor, constructible from a pointer to the derived class and 
 *   that behaves like a distance functor from the begin() iterator to a given point.
 *   (see SimpleRandomAccessRangeFromPoint and SimpleRandomAccessConstRangeFromPoint)
 */
template < typename TDerived >
class IteratorCompletionTraits;

/// Class that uses CRTP to add reverse iterators and ranges to a derived class.
template < 
  typename TDerived
>
class IteratorCompletion
  {
  public:

    using Iterator        = typename IteratorCompletionTraits<TDerived>::Iterator;
    using ConstIterator   = typename IteratorCompletionTraits<TDerived>::ConstIterator;
    using DistanceFunctor = typename IteratorCompletionTraits<TDerived>::DistanceFunctor;

    using ReverseIterator       = boost::reverse_iterator<Iterator>;
    using ConstReverseIterator  = boost::reverse_iterator<ConstIterator>;
    using Range       = SimpleRandomAccessRangeFromPoint< ConstIterator, Iterator, DistanceFunctor >;
    using ConstRange  = SimpleRandomAccessConstRangeFromPoint< ConstIterator, DistanceFunctor >;
    using Difference  = std::ptrdiff_t;

    /**
     * @return  a mutable reverse-iterator pointing to the last value.
     * @warning the derived class must have a begin() method.
     */
    ReverseIterator rbegin()
      {
        return ReverseIterator{ static_cast<TDerived*>(this)->end() };
      }

    /**
     * @return  a constant reverse-iterator pointing to the last value.
     * @warning the derived class must have a cend() method.
     */
    inline
    ConstReverseIterator rbegin() const
      {
        return ConstReverseIterator{ static_cast<TDerived*>(this)->end() };
      }
    
    /**
     * @return  a constant reverse-iterator pointing to the last value (C++11).
     * @warning the derived class must have a cend() method.
     */
    inline
    ConstReverseIterator crbegin() const
      {
        return ConstReverseIterator{ static_cast<TDerived*>(this)->end() };
      }
    
    /**
     * @return  a mutable reverse-iterator pointing before the first value.
     * @warning the derived class must have a end() method.
     */
    inline
    ReverseIterator rend()
      {
        return ReverseIterator{ static_cast<TDerived*>(this)->begin() };
      }

    /**
     * @return  a constant reverse-iterator pointing before the first value.
     * @warning the derived class must have a cbegin() method.
     */
    inline
    ConstReverseIterator rend() const
      {
        return ConstReverseIterator{ static_cast<TDerived*>(this)->begin() };
      }
    
    /**
     * @return  a constant reverse-iterator pointing before the first value (C++11).
     * @warning the derived class must have a cbegin() method.
     */
    inline
    ConstReverseIterator crend() const
      {
        return ConstReverseIterator{ static_cast<TDerived*>(this)->begin() };
      }

    /**
     * @return  a mutable range over the derived class values.
     * @warning the derived class must have mutable iterators and a distance functor.
     */
    inline
    Range range()
      {
        TDerived* const derived = static_cast<TDerived*>(this);
        return { 
            derived->begin(), 
            derived->end(),
            typename IteratorCompletionTraits<TDerived>::DistanceFunctor( derived )
        };
      }

    /**
     * @return  a constant range over the derived class values.
     * @warning the derived class must have constant iterators and a distance functor.
     */
    inline
    ConstRange constRange() const
      {
        TDerived const* const derived = static_cast<TDerived const*>(this);
        return { 
            derived->begin(), 
            derived->end(), 
            typename IteratorCompletionTraits<TDerived>::DistanceFunctor( derived )
        };
      }

  protected:

    /// Protected destructor to avoid memory leak.
    ~IteratorCompletion()
      {}
  };

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

