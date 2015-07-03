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
 * @file MultiPhaseField2.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @author Roland Denis (\c roland.denis@univ-smb.fr )
 * LAboratory of MAthematics - LAMA (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/04/28
 *
 * This file is part of the DGtal library.
 */

#if defined(MultiPhaseField2_RECURSES)
#error Recursive header files inclusion detected in MultiPhaseField2.h
#else // defined(MultiPhaseField2_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MultiPhaseField2_RECURSES

#if !defined MultiPhaseField2_h
/** Prevents repeated inclusion of headers. */
#define MultiPhaseField2_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cstddef>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"

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
  // template class MultiPhaseField2
  /**
   * Description of template class 'MultiPhaseField2' <p>
   * \brief Aim: This class is a way of deforming an image of labels. 
   * Each region (ie. set of points having a same label) is viewed as 
   * the set of points having a value greater than 0.5 for a given phase field. 
   * Each region is evolved through its phase field. 
   *
   * @tparam TLabelImage a model of CImage (storing labels)
   * @tparam TFieldImage a model of CImage (storing phase field values)
   * @tparam TEvolver a model of phase field evolver
   */
  template <typename TLabelImage, typename TFieldImage, typename TMultiFieldImage>
  class MultiPhaseField2;

  template <
    typename TLabelImage,
    typename TFieldImage,
    typename TDomain, typename TContainer, typename TApproximation, typename TBoundingBox
  >
  class MultiPhaseField2< TLabelImage, TFieldImage, ApproximatedMultiImage<TDomain, TContainer, TApproximation, TBoundingBox> >
  {

    // ----------------------- Types check -----------------------

    BOOST_CONCEPT_ASSERT(( concepts::CImage<TLabelImage> )); 
    BOOST_CONCEPT_ASSERT(( concepts::CImage<TFieldImage> )); 
    BOOST_STATIC_ASSERT
    (( concepts::ConceptUtils::SameType< typename TLabelImage::Point,
       typename TFieldImage::Point>::value ));
    BOOST_STATIC_ASSERT(( concepts::ConceptUtils::SameType< typename TFieldImage::Domain, TDomain >::value ));


    // ----------------------- Types ------------------------------
  public:

    /// Image of labels
    typedef TLabelImage LabelImage;
    typedef typename LabelImage::Value Label;
    typedef typename LabelImage::Domain Domain;
    typedef typename Domain::Point Point;

    /// Images of phase field values
    typedef TFieldImage FieldImage;
    typedef typename TFieldImage::Value Value;

    using MultiImage = ApproximatedMultiImage<TDomain, TContainer, TApproximation, TBoundingBox>;

    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:


    // ------------------------- References --------------------------------
    /**
     * Reference on the image of labels
     */
    LabelImage& myLabelImage; 


    // ------------------------- Data --------------------------------
   
    /** Container of the fields
     */
    MultiImage myFields;
    
    /**
     * List of labels
     */
    std::vector<Label> myLabels; 

    //! Initial volumes
    std::vector<Value> myInitVolume;

    //! Epsilon 
    Value myEpsilon; 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aI an image of labels
     */
    MultiPhaseField2(LabelImage& aI, Value epsilon, bool calcDistance = true);

    /**
     * Destructor. Does nothing.
     */
    ~MultiPhaseField2();


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
    std::size_t updateLabels();

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
     * Return the number of phases.
     */
    size_t getNumPhase() const; 

    /**
     * Return one of the phase fields
     */
    FieldImage getPhase( size_t i ) const;

    /**
     * Return the perimeter of each phase fields
     */
    std::vector<Value> getPerimeters() const;

    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    MultiPhaseField2 ( const MultiPhaseField2 & other ) = delete;

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    MultiPhaseField2 & operator= ( const MultiPhaseField2 & other ) = delete;

  private:



    // ------------------------- Internals ------------------------------------
  private:

    /**
     * Init @a aImage by a signed distance to the frontier
     * of the region having @a aLabel as label
     * @param aLabel region id
     * @param aImage image to initialize
     */
    void getSignedDistance(const Label& aLabel, FieldImage& aImage) const; 

  }; // end of class MultiPhaseField2


  /**
   * Overloads 'operator<<' for displaying objects of class 'MultiPhaseField2'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'MultiPhaseField2' to write.
   * @return the output stream after the writing.
   */
  template <
    typename TLabelImage,
    typename TFieldImage,
    typename TDomain, typename TContainer, typename TApproximation, typename TBoundingBox
  >
  std::ostream&
  operator<< ( std::ostream & out, 
		    DGtal::MultiPhaseField2< TLabelImage, TFieldImage, DGtal::ApproximatedMultiImage<TDomain, TContainer, TApproximation, TBoundingBox> > const& object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "MultiPhaseField2.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MultiPhaseField2_h

#undef MultiPhaseField2_RECURSES
#endif // else defined(MultiPhaseField2_RECURSES)

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

