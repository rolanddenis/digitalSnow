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
 * @file NoBoundingBox.h
 * @author Roland Denis (\c roland.denis@univ-savoie.fr )
 * LAboratory of MAthematics - LAMA (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/03/26
 *
 * This file is part of the DGtal library.
 */

#if defined(NoBoundingBox_RECURSES)
#error Recursive header files inclusion detected in NoBoundingBox.h
#else // defined(NoBoundingBox_RECURSES)
/** Prevents recursive inclusion of headers. */
#define NoBoundingBox_RECURSES

#if !defined NoBoundingBox_h
/** Prevents repeated inclusion of headers. */
#define NoBoundingBox_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class NoBoundingBox
  /**
   * @brief Implements no Bounding Box.
   *
   * @tparam TDomain  a HyperRectDomain.
   */
  template <
    typename TDomain
  >
  class NoBoundingBox;

  // Implementation for HyperRectDomain
  template <
    typename TSpace
  >
  class NoBoundingBox< HyperRectDomain<TSpace> >
  {
    // ----------------------- Typedefs & constants ---------------------------
  public:
    typedef HyperRectDomain<TSpace> Domain;
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
    NoBoundingBox( Domain const& domain ) : m_domain(domain) {}
    
    /**
     * Destructor.
     */
    ~NoBoundingBox() {}

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Reset point counters and bounding box.
     */
    static inline void reset () {}

    /**
     * Add a point in the bounding box.
     * @param point   the point to add.
     */
    static inline void addPoint ( Point const& /* point */ ) {}

    /**
     * Remove a point from the bounding box.
     * @param point   the point to remove.
     */
    static inline void removePoint ( Point const& /* point */ ) {}

    /**
     * Get lower bound of the bounding box.
     */
    inline Point lowerBound () const { return m_domain.lowerBound(); }

    /**
     * Get upper bound of the bounding box.
     */
    inline Point upperBound () const { return m_domain.upperBound(); }

    /**
     * Get bounding box as a HyperRectDomain with buffer zone.
     * @param buffer The buffer size as a point.
     */
    inline Domain 
    getBoundingBox( Point const& buffer = Point::diagonal(0) ) const
      {
        return Domain{ lowerBound() - buffer, upperBound() + buffer };
      }

    /**
     * Returns true if the bounding box is empty
     */
    static inline bool isEmpty() { return false; } ///< \todo check if domain is empty

    /**
     * Returns true if the bounding box is not empty
     */
    static inline bool isNotEmpty() { return true; }
    
    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    inline void selfDisplay ( std::ostream & out ) const
      {
        m_domain.selfDisplay( out );
      }

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    static inline bool isValid() { return true; }

    // ------------------------- Private Datas --------------------------------
  private:
    Domain m_domain;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    NoBoundingBox();

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class NoBoundingBox


  /**
   * Overloads 'operator<<' for displaying objects of class 'NoBoundingBox'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'NoBoundingBox' to write.
   * @return the output stream after the writing.
   */
  template <typename T>
  std::ostream&
  operator<< ( std::ostream & out, const NoBoundingBox<T> & object )
    {
      object.selfDisplay(out);
      return out;
    }

} // namespace DGtal


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NoBoundingBox_h

#undef NoBoundingBox_RECURSES
#endif // else defined(NoBoundingBox_RECURSES)

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

