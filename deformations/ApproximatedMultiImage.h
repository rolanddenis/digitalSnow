#pragma once

#include <boost/concept_check.hpp>
#include <cstddef>
#include <vector>
#include <limits>
#include <ostream>
#include <cmath>
#include <utility>

#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/domains/CDomain.h>
#include "ValueApproximations.h"
#include "NoBoundingBox.h"
#include "ImageView.h"
#include "Linearizer.h"

namespace approximated_multi_image
{

  /**
   * Type traits to get the value type of a container.
   *
   * For most containers, the typedef value_type is what we want but others containers
   * like std::map of DGtal::LabelledMap exhibits instead the mapped_type typedef 
   * ( value_type represents, in this cases, the pair key|value )
   *
   * @tparam TContainer Type of the container
   */
  template < typename TContainer >
  struct ValueType
    {
      typedef typename TContainer::value_type type; ///< The type of the contained data.
    };

  /**
   * Specialization of ValueType for LabelledMap.
   */
  template < typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M >
  struct ValueType< DGtal::LabelledMap<TData, L, TWord, N, M> >
    {
      typedef TData type;
    };


} // namespace approximated_multi_image

namespace DGtal
{

/// \todo move linearizer from ImageContainerBySTLVector into anonymous namespace to avoid polutting DGtal namespace.
using namespace approximated_multi_image;

/** Informations about multiple image container
 */
template <typename t>
struct MultiImageInfos;

/// Memory usage of multiple LabelledMap based on the size histogram.
template <
  typename TData, typename TWord,
  typename TDomain,
  typename TSizeHist
>
size_t getMultiImageMemoryUsage(
    TDomain const& aDomain,
    TSizeHist const& aSizeHist,
    unsigned int L, unsigned int N, unsigned int M
);

/// Find best parameters for LabelledMap given the size histogram.
template <
  typename TData, typename TWord,
  typename TDomain,
  typename TSizeHist
>
std::pair<unsigned int, unsigned int> getOptimalLabelledMap( TDomain const& aDomain, TSizeHist const& aSizeHist, unsigned int L );

/**
 * Multiple images container with approximation and bounding box capabilities.
 *
 * @tparam TDomain    The domain of the images.
 * @tparam TContainer The container used to stored, at each point of the domain, the values of the images.
 * @tparam TApproximation The predicate used to approximate the values.
 * @tparam TBoundingBox   The type of bounding box for the non-approximated values.
 */
template <
  typename TDomain,
  typename TContainer,
  typename TApproximation = approximations::NoValueApproximation< typename ValueType<TContainer>::type >,
  typename TBoundingBox = NoBoundingBox<TDomain>
>
class ApproximatedMultiImage;

/**
 * Specialization of ApproximatedMultiImage for LabelledMap container.
 */
template <
  typename TDomain,
  typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M,
  typename TApproximation,
  typename TBoundingBox
>
class ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox >
  {
  public:

    BOOST_CONCEPT_ASSERT( (concepts::CDomain<TDomain>) );
    BOOST_CONCEPT_ASSERT( (concepts::CValueApproximation<TApproximation>) );
    static_assert( std::is_same< typename TApproximation::value_type, TData>::value, "Container and approximation must have same value_type");
    /// \todo Concept check for bounding box

    // Typedefs
    class Reference;

    // DGtal typedefs
    using Self          = ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox >;
    using Domain        = TDomain; ///< Type of the domain.
    using Container     = DGtal::LabelledMap<TData, L, TWord, N, M>; ///< Type of the values container.
    using Approximation = TApproximation; ///< Type of the approximation.
    using BoundingBox   = TBoundingBox;   ///< Type of the bounding box.
    using Label         = typename Container::Label;  ///< Type of the label used to identify an image.
    using Value         = TData;  ///< Type of the stored values.
    using Point         = typename TDomain::Point;  ///< Type of a point in the domain.
    using Size          = typename TDomain::Size;   ///< Type of the linearized index of a point.
    using ValueConstIterator = typename Container::ConstIterator;   ///< Type of the const-iterator over the values associated to a given point.
    using FullImage     = ImageView< Self, image_view::FullDomain >; ///< Full image view.
    using ConstFullImage = ImageView< const Self, image_view::FullDomain >; ///< Full image constant view;
    using BBImage       = ImageView< Self, image_view::BoundingBoxAsDomain >; ///< Image view from bounding box.
    using ConstBBImage  = ImageView< const Self, image_view::BoundingBoxAsDomain >; ///< Image constant view from bounding box.


    // STL typedefs
    // ...

    // Friend with image view
    friend FullImage;
    friend ConstFullImage;
    friend BBImage;
    friend ConstBBImage;

    /**
     * Constructor.
     */
    ApproximatedMultiImage( Domain const& aDomain, Approximation const& anApprox = Approximation{} ) noexcept
      : myDomain{aDomain}, myImages{aDomain.size()}, 
        myApproximation{anApprox}, myBoundingBoxes{L, BoundingBox{aDomain}},
        myExtent{ aDomain.upperBound() - aDomain.lowerBound() + Point::diagonal(1) }
    {}

    /**
     * Get a value given the linearized index of the point.
     *
     * @param anIndex The index of the point.
     * @param aLabel  The label of the image.
     */
    Value getValueByIndex( Size anIndex, Label aLabel ) const noexcept
      {
        Container const& values = myImages[ anIndex ];
        
        if ( values.count(aLabel) > 0 )
          {
            return values.fastAt(aLabel);
          } 
        else
          {
            return myApproximation.default_value;
          }
      }

    /**
     * Get a value.
     *
     * @param aPoint The point of the domain.
     * @param aLabel The label of the image from which to get the value.
     */
    inline
    Value getValue( Point const& aPoint, Label aLabel ) const noexcept
      {
        return getValueByIndex( linearized(aPoint), aLabel );
      }


        
    /**
     * Set a value.
     *
     * @param aPoint  The point of the domain.
     * @param aLabel  The label of the image in which to set the value.
     * @param aValue  The value to be set.
     * @param anIndex The linearized index of the point (by default, calculated from the given point).
     */
    inline
    void setValue( Point const& aPoint, Label aLabel, Value aValue ) noexcept
      {
        setValue( aPoint, aLabel, aValue, linearized(aPoint) );
      }

    inline
    void setValue( Point const& aPoint, Label aLabel, Value aValue, Size anIndex ) noexcept
      {
        Container & values = myImages[ anIndex ];

        if ( myApproximation.eval( aValue ) )
          {
            const auto it = values.find(aLabel);
            if ( it != values.end() )
              {
                values.erase(it);
                myBoundingBoxes[aLabel].removePoint(aPoint);
              }
          }
        else
          {
            if ( values.count(aLabel) > 0 )
              {
                values.fastAt(aLabel) = aValue;
              }
            else
              {
                values[aLabel] = aValue;
                myBoundingBoxes[aLabel].addPoint(aPoint);
              }
          }
      }

    /**
     * Return the domain of the images.
     */
    inline
    Domain const& domain() const noexcept
      {
        return myDomain;
      }
    
    /**
     * Return the bounding box (with buffer) of an image.
     * @param aLabel  Label of the image.
     * @param buffer  Buffer around the initial bounding box.
     */
    inline
    Domain getBoundingBox( Label aLabel, Point const& buffer = Point::diagonal(0) ) const noexcept
      {
        return myBoundingBoxes[ aLabel ].getBoundingBox( buffer );
      }
    
    /**
     * Return the value container associated to a point.
     * \todo mutable version.
     *
     * @param aPoint The point.
     */
    inline
    Container const& operator() ( Point const& aPoint ) const noexcept
      {
        return myImages[ linearized(aPoint) ];
      }

    /**
     * Return a full-domain view of the image with the given label
     * @param aLabel The image label.
     */
    inline
    FullImage operator[] ( Label aLabel ) noexcept
      {
        return { *this, aLabel };
      }

    /**
     * Return a full-domain constant view of the image with the given label
     * @param aLabel The image label.
     */
    inline
    ConstFullImage operator[] ( Label aLabel ) const noexcept
      {
        return { *this, aLabel };
      }

    /**
     * Return an image view restricted to his bounding box (with buffer).
     * @param aLabel  The image label.
     * @param buffer  The buffer around the bounding box.
     */
    inline
    BBImage getBBImage( Label aLabel, Point buffer = Point::diagonal(0) )
      {
        BBImage image{ *this, aLabel };
        image.buffer() = buffer;
        return image;
      }
    
    /**
     * Return an image constant-view restricted to his bounding box (with buffer).
     * @param aLabel  The image label.
     * @param buffer  The buffer around the bounding box.
     */
    inline
    ConstBBImage getBBImage( Label aLabel, Point buffer = Point::diagonal(0) ) const
      {
        ConstBBImage image{ *this, aLabel };
        image.buffer() = buffer;
        return image;
      }


    /// Get informations about this class
    MultiImageInfos<Self> getInfos() const
      {
        using std::size_t;
        MultiImageInfos<Self> infos;
        infos.labelMin = std::numeric_limits<size_t>::max();
        infos.labelMax = std::numeric_limits<size_t>::min();
        size_t label_sum = 0;
        size_t label_sqr_sum = 0;
        size_t cnt = 0;
        std::array<size_t, L> support; support.fill(0);
        infos.imageVolume.fill(0);
        infos.labelHist.fill(0);

        for ( auto const& point_values : myImages )
          {
            const size_t size = point_values.size();
            infos.labelMin = std::min( infos.labelMin, size );
            infos.labelMax = std::max( infos.labelMax, size );
            label_sum += size;
            label_sqr_sum += size*size;
            ++infos.labelHist[size];

            for ( auto const& label_value : point_values )
              {
                ++(support[label_value.first]);
                infos.imageVolume[label_value.first] += label_value.second;
              }

            ++cnt;
          }

        infos.labelMean = double(label_sum) / cnt;
        infos.labelSDeviation = std::sqrt( double(label_sqr_sum)/cnt - std::pow(infos.labelMean, 2) );

        for ( size_t i = 0; i < L; ++i )
          {
            infos.imageSupport[i] = double(support[i]) / cnt;
            infos.imageBB[i] = double( getBoundingBox(i).size() ) / cnt;
          }

        infos.memoryUsage = getMultiImageMemoryUsage<TData,TWord>( myDomain, infos.labelHist, L, N, M );
        auto const bestSettings = getOptimalLabelledMap<TData,TWord>( myDomain, infos.labelHist, L );
        infos.bestN = bestSettings.first;
        infos.bestM = bestSettings.second;
        infos.bestMemoryUsage = getMultiImageMemoryUsage<TData,TWord>( myDomain, infos.labelHist, L, infos.bestN, infos.bestM );

        return infos;
      }

    /// Reference to an approximated value.
    class Reference
      {
      public:
        Reference( Self& aMultiImage, Point aPoint, Label aLabel, Size anIndex )
          : myMultiImage{ aMultiImage }, myPoint{ aPoint }, myLabel{ aLabel }, myIndex{ anIndex }
        {}

        Reference( Self& aMultiImage, Point aPoint, Label aLabel )
          : Reference{ aMultiImage, aPoint, aLabel, aMultiImage.linearized(aPoint) }
        {}

        operator Value() const
          {
            return myMultiImage.getValueByIndex( myIndex, myLabel );
          }

        Reference& operator= ( Value aValue )
          {
            myMultiImage.setValue( myPoint, myLabel, aValue, myIndex );
            return *this;
          }

      private:
        Self& myMultiImage;
        Point myPoint;
        Label myLabel;
        Size myIndex;
      };

  public:

    /**
     * Return the linearized index of a point.
     *
     * @param aPoint The point.
     */
    inline
    Size linearized( Point const& aPoint ) const
      {
        return Linearizer<Domain,ColMajorStorage>::apply(aPoint, myDomain.lowerBound(), myExtent);
      }


  private:
    Domain myDomain;
    std::vector<Container> myImages;
    Approximation myApproximation;
    std::vector<BoundingBox> myBoundingBoxes;
    Point myExtent;
    
  }; // class ApproximatedMultiImage

  /// Informations about this class
  template <
    typename TDomain,
    typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M,
    typename TApproximation,
    typename TBoundingBox
  >
  struct MultiImageInfos< DGtal::ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox > > 
    {
      std::size_t labelMin, labelMax; // Minimum and maximum number of labels stored in each point.
      double labelMean, labelSDeviation; // Mean and standard deviation of the number of labels.
      std::array<size_t, L+1> labelHist; // Histogram of the number of labels stored.
      std::array<TData, L> imageVolume; // Volume of each image.
      std::array<double, L> imageSupport; // Relative support of each image (in [0,1]).
      std::array<double, L> imageBB; // Relative volume of the bounding box of each image (in [0,1]).
      size_t memoryUsage; // Memory usage of this structure
      unsigned int bestN, bestM; // Optimal LabelledMap settings.
      size_t bestMemoryUsage; // Memory usage for the optimal LabelledMap settings.
    };

  /// Memory usage of multiple LabelledMap based on the size histogram.
  template <
    typename TData, typename TWord,
    typename TDomain,
    typename TSizeHist
  >
  size_t getMultiImageMemoryUsage(
      TDomain const& aDomain,
      TSizeHist const& aSizeHist,
      unsigned int L, unsigned int N, unsigned int M
  )
    {
      const size_t sizeof_LabelledMap = 
            sizeof(TWord) * ( L/(8*sizeof(TWord)) + ( L % 8*sizeof(TWord) == 0 ? 0 : 1 ) )
          + N * sizeof(TData) + std::max(sizeof(TData), sizeof(TData*));

      const size_t sizeof_LastBlock = M * sizeof(TData) + sizeof(TData*);
      
      size_t usage = 0;
      for ( unsigned int i = N+2; i <= L; ++i )
        usage += aSizeHist[i] * ( (i-N)/M + ( (i-N) % M == 0 ? 0 : 1 ) );


      usage = sizeof_LastBlock * usage + aDomain.size() * sizeof_LabelledMap;

      return usage;
    }

  /// Find best parameters for LabelledMap given the size histogram.
  template <
    typename TData, typename TWord,
    typename TDomain,
    typename TSizeHist
  >
  std::pair<unsigned int, unsigned int> 
      getOptimalLabelledMap( TDomain const& aDomain, TSizeHist const& aSizeHist, unsigned int L )
    {
      std::pair<unsigned int, unsigned int> param = { 1, 1 };
      const size_t domain_size = aDomain.size();
      size_t min_mem = getMultiImageMemoryUsage<TData,TWord>(aDomain, aSizeHist, L, param.first, param.second);
      size_t max_N = min_mem/sizeof(TData) + ( min_mem % sizeof(TData) == 0 ? 0 : 1 );

      for ( unsigned int N = 1; N < L && N < max_N; ++N )
        {
          size_t last_mem = std::numeric_limits<size_t>::max();
          for ( unsigned int M = 1; M < L; ++M )
            {
              const size_t mem = getMultiImageMemoryUsage<TData,TWord>(aDomain, aSizeHist, L, N, M);
              if ( mem >= last_mem ) break;
              last_mem = mem;
              if ( mem < min_mem || ( mem == min_mem && (N > param.first || M > param.second) ) )
                {
                  min_mem = mem;
                  param.first = N;
                  param.second = M;
                }
            }
          max_N = min_mem/sizeof(TData) + ( min_mem % sizeof(TData) == 0 ? 0 : 1 );
        }

      return param;
    }

  /// Display of ApproximatedMultiImage informations
  //
  template <
    typename TDomain,
    typename TData, unsigned int L, typename TWord, unsigned int N, unsigned int M,
    typename TApproximation,
    typename TBoundingBox
  >
  std::ostream& operator<< ( std::ostream& out, typename DGtal::MultiImageInfos< DGtal::ApproximatedMultiImage< TDomain, DGtal::LabelledMap<TData, L, TWord, N, M>, TApproximation, TBoundingBox > > const& infos )
    {
      out << "Label count: "
          << "min=" << infos.labelMin << " ; "
          << "max=" << infos.labelMax << " ; "
          << "mean=" << infos.labelMean << " ; "
          << "sdev=" << infos.labelSDeviation
          << std::endl;

      out << "Best settings: N=" << infos.bestN << "(" << N 
          << ") ; M=" << infos.bestM << "(" << M 
          << ") ; memory=" << infos.bestMemoryUsage << "(" << infos.memoryUsage << ")"
          << std::endl;

      bool first_label = true;
      for ( std::size_t i = 0; i < L; ++i ) 
        {
          if ( infos.imageVolume[i] > 0 )
            {
              if (! first_label) out << " ; ";
              out << "#" << i
                  << " V"  << infos.imageVolume[i]
                  << " S"  << std::fixed << std::setprecision(1) << 100*infos.imageSupport[i]
                  << " BB" << 100*infos.imageBB[i]
                  << " R"  << std::setprecision(2) << infos.imageBB[i]/infos.imageSupport[i];
              first_label = false;
            }
        }
      out << std::endl;

      return out;
    }

} // namespace DGtal


/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

