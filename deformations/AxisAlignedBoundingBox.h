/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file AxisAlignedBoundingBox.h
 * @author Roland Denis (\c roland.denis@univ-savoie.fr )
 * LAboratory of MAthematics - LAMA (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/03/20
 *
 * This file is part of the DGtal library.
 */

#if defined(AxisAlignedBoundingBox_RECURSES)
#error Recursive header files inclusion detected in AxisAlignedBoundingBox.h
#else // defined(AxisAlignedBoundingBox_RECURSES)
/** Prevents recursive inclusion of headers. */
#define AxisAlignedBoundingBox_RECURSES

#if !defined AxisAlignedBoundingBox_h
/** Prevents repeated inclusion of headers. */
#define AxisAlignedBoundingBox_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <array>
#include <vector>

#include <DGtal/base/Common.h>
#include <DGtal/kernel/CIntegralNumber.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class AxisAlignedBoundingBox
  /**
   * @brief Axis Aligned Bounding Box on HyperRectDomain
   *
   * @tparam TDomain    a HyperRectDomain.
   * @tparam TCounter   type of the pointer counter in each dimension.
   * @tparam N...       the dimensions that must be taken into account. All dimensions by default.
   */
  template <
    typename TDomain,
    typename TCounter,
    DGtal::Dimension ...N
  >
  class AxisAlignedBoundingBox;

  // Fallback for HyperRectDomain if no working dimensions are specified
  template <
    typename TSpace,
    typename TCounter
  >
  class AxisAlignedBoundingBox< HyperRectDomain<TSpace>, TCounter >;
  
  // Implementation for HyperRectDomain
  template <
    typename TSpace,
    typename TCounter,
    DGtal::Dimension ...N
  >
  class AxisAlignedBoundingBox< HyperRectDomain<TSpace>, TCounter, N... >
  {
    // ----------------------- Concepts check ---------------------------------
    BOOST_CONCEPT_ASSERT(( concepts::CIntegralNumber<TCounter> ));
    /// \todo check N...

    // ----------------------- Typedefs & constants ---------------------------
  public:
    typedef HyperRectDomain<TSpace> Domain;
    typedef TCounter          Counter;
    typedef typename Domain::Point     Point;
    typedef typename Domain::Dimension Dimension;

    BOOST_STATIC_CONSTANT( Dimension, dimension = Domain::dimension );

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param domain  domain containing the bounding box.
     */
    AxisAlignedBoundingBox( Domain const& domain );

    /**
     * Destructor.
     */
    ~AxisAlignedBoundingBox();

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Reset point counters and bounding box.
     */
    void reset ();

    /**
     * Add a point in the bounding box.
     * @param point   the point to add.
     */
    void addPoint ( Point const& point );

    /**
     * Remove a point from the bounding box.
     * @param point   the point to remove.
     */
    void removePoint ( Point const& point );

    /**
     * Get lower bound of the bounding box.
     */
    Point lowerBound () const;

    /**
     * Get upper bound of the bounding box.
     */
    Point upperBound () const;

    /**
     * Get bounding box as a HyperRectDomain.
     * @param 
     */
    Domain getBoundingBox( Point const& buffer = Point::diagonal(0) ) const;

    /**
     * Returns true if the bounding box is empty
     */
    bool isEmpty() const;

    /**
     * Returns true if the bounding box is not empty
     */
    bool isNotEmpty() const;
    
    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    static bool isValid();

    // ------------------------- Private Datas --------------------------------
  private:
    Domain m_domain; ///< \remark We do not really need the domain (that take 88bytes!) ... only his lower & upper bounds. 
    std::array< std::vector<Counter>, sizeof...(N) > m_counters;
    Point m_lowerBound, m_upperBound;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    AxisAlignedBoundingBox();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    AxisAlignedBoundingBox ( const AxisAlignedBoundingBox & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    AxisAlignedBoundingBox & operator= ( const AxisAlignedBoundingBox & other );

    // ------------------------- Internals ------------------------------------
  private:

    template < std::size_t I, DGtal::Dimension ...M > struct reset_impl;
    template < std::size_t I, DGtal::Dimension ...M > struct addPoint_impl;
    template < std::size_t I, DGtal::Dimension ...M > struct removePoint_impl;
    template < std::size_t I, DGtal::Dimension ...M > struct isEmpty_impl;

  }; // end of class AxisAlignedBoundingBox


  /**
   * Overloads 'operator<<' for displaying objects of class 'AxisAlignedBoundingBox'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'AxisAlignedBoundingBox' to write.
   * @return the output stream after the writing.
   */
  template < typename TDomain, typename TCounter, DGtal::Dimension ...N >
  std::ostream&
  operator<< ( std::ostream & out, const AxisAlignedBoundingBox<TDomain, TCounter, N...> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "AxisAlignedBoundingBox.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined AxisAlignedBoundingBox_h

#undef AxisAlignedBoundingBox_RECURSES
#endif // else defined(AxisAlignedBoundingBox_RECURSES)

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

