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
 * @file MultiPhaseFieldDrip.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @author Roland Denis (\c roland.denis@univ-smb.fr )
 * LAboratory of MAthematics - LAMA (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/04/28
 *
 * This file is part of the DGtal library.
 */

#if defined(MultiPhaseFieldDrip_RECURSES)
#error Recursive header files inclusion detected in MultiPhaseFieldDrip.h
#else // defined(MultiPhaseFieldDrip_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MultiPhaseFieldDrip_RECURSES

#if !defined MultiPhaseFieldDrip_h
/** Prevents repeated inclusion of headers. */
#define MultiPhaseFieldDrip_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cstddef>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"
#include "DGtal/kernel/PointVector.h"

#include "DGtal/base/CowPtr.h"
#include "DGtal/base/CountedPtr.h"


//for getSignedDistance private method
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/geometry/volumes/distance/DistanceTransformation.h>

#include "deformationFunctions.h"
#include "ApproximatedMultiImage.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{


  /////////////////////////////////////////////////////////////////////////////
  // template class MultiPhaseFieldDrip
  /**
   * Description of template class 'MultiPhaseFieldDrip' <p>
   * \brief Aim: This class is a way of deforming an image of labels.
   * Each region (ie. set of points having a same label) is viewed as
   * the set of points having a value greater than 0.5 for a given phase field.
   * Each region is evolved through its phase field.
   *
   * @tparam TFieldImage a model of CImage (storing phase field values)
   * @tparam TEvolver a model of phase field evolver
   */
  template < typename TFieldImage, typename TMultiFieldImage>
  class MultiPhaseFieldDrip;

  template <
    typename TFieldImage,
    typename TDomain, typename TContainer, typename TApproximation, typename TBoundingBox
  >
  class MultiPhaseFieldDrip< TFieldImage, ApproximatedMultiImage<TDomain, TContainer, TApproximation, TBoundingBox> >
  {

    // ----------------------- Types check -----------------------

    BOOST_CONCEPT_ASSERT(( concepts::CImage<TFieldImage> ));
    BOOST_STATIC_ASSERT(( concepts::ConceptUtils::SameType< typename TFieldImage::Domain, TDomain >::value ));


    // ----------------------- Types ------------------------------
  public:

    /// Images of phase field values
    using FieldImage  = TFieldImage;
    using Value       = typename FieldImage::Value;
    using Domain      = typename FieldImage::Domain;
    using Point       = typename Domain::Point;
    using RealPoint   = typename TDomain::Space::RealPoint;
    using MultiImage  = ApproximatedMultiImage<TDomain, TContainer, TApproximation, TBoundingBox>;

    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Data --------------------------------

    /// Max number of phases.
    std::size_t myMaxPhaseCnt;

    /** Container of the fields
     */
    MultiImage myFields;

    //! Initial volumes
    std::vector<Value> myTargetVolume;

    //! Epsilon
    Value myEpsilon;

  public:
    RealPoint myRealExtent = RealPoint::diagonal(1);

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aI an image of labels
     */
    MultiPhaseFieldDrip( Domain const& aDomain, std::size_t maxPhaseCnt, Value epsilon);

    /**
     * Destructor. Does nothing.
     */
    ~MultiPhaseFieldDrip();


    /**
     * Deform the image of labels during @a aT
     *
     * @param aT time step
     * @return time spent during the deformation
     * (equal to aT).
     */
    double update(const double& aT);

    /**
     * Updates image of labels and returns modification count.
     */
    template < typename TLabelImage >
    std::size_t updateLabels( TLabelImage & aLabelImage );

    /**
     * Calculates and diplays informations about the phase fields.
     */
    void dispInfos() const;


    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;


    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Return the current number of phases.
     */
    size_t getNumPhase() const;

    /**
     * Return the maximum number of phases.
     */
    size_t getMaxNumPhase() const;

    /**
     * Return one of the phase fields
     */
    FieldImage getPhase( size_t i ) const;

    /** Adds a new phase in the domain is the maximum number is not reached
     * @return true if a new phase has been added.
     */
    bool addPhase();

    template < typename TGenerator >
    bool addPhase( TGenerator & gen );

    /**
     * Return the perimeter of each phase fields
     */
    std::vector<Value> getPerimeters() const;

    /// Returns a constant reference on the fields container.
    MultiImage const& getPhasesContainer() const
      {
        return myFields;
      }

    void updateDomainSize();

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    MultiPhaseFieldDrip ( const MultiPhaseFieldDrip & other ) = delete;

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    MultiPhaseFieldDrip & operator= ( const MultiPhaseFieldDrip & other ) = delete;

  }; // end of class MultiPhaseFieldDrip


  /**
   * Overloads 'operator<<' for displaying objects of class 'MultiPhaseFieldDrip'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'MultiPhaseFieldDrip' to write.
   * @return the output stream after the writing.
   */
  template <
    typename TFieldImage,
    typename TDomain, typename TContainer, typename TApproximation, typename TBoundingBox
  >
  std::ostream&
  operator<< ( std::ostream & out,
      DGtal::MultiPhaseFieldDrip< TFieldImage, DGtal::ApproximatedMultiImage<TDomain, TContainer, TApproximation, TBoundingBox> > const& object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "MultiPhaseFieldDrip.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MultiPhaseFieldDrip_h

#undef MultiPhaseFieldDrip_RECURSES
#endif // else defined(MultiPhaseFieldDrip_RECURSES)

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
