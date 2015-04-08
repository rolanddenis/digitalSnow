#include <iostream>
#include <cstddef>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#include <DGtal/base/Common.h>
#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include <boost/math/constants/constants.hpp>

#include "AxisAlignedBoundingBox.h"
#include "ValueApproximations.h"
#include "ApproximatedMultiImage.h"
#include "AxisAlignedBoundingBox.h"
#include "ImageView.h"

using std::size_t;
using namespace DGtal;

constexpr 
size_t ipow( size_t base, size_t exp )
{
  return exp > 0 ? base * ipow(base, exp-1) : 1;
}

template <
  typename TImage,
  size_t N
>
class BenchVectorOfImages
{
public:

  using Image = TImage;
  using Domain = typename Image::Domain;
  using Value = typename Image::Value;
  using Point = typename Domain::Point;

  BenchVectorOfImages( Domain const& domain )
    {
      const size_t L = ipow( N, Domain::dimension );
      myImages.reserve(L);
      for ( size_t i = 0; i < L; ++i )
        myImages.emplace_back(domain);
    }
  
  inline size_t nImages() const
    {
      return ipow(N, Domain::dimension);
    }

  Value generateData( Value radius, Value eps )
    {
      const size_t dimension = Domain::dimension;
      Domain const& domain = myImages[0].domain();

      // Shift between each phase center
      const Point extent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
      Value dx[dimension];
      for (size_t i = 0; i < dimension; ++i)
        dx[i] = Value(extent[i]) / N;

      Domain pos_domain{ Point::diagonal(0), Point::diagonal(N-1) };

      size_t Label = 0;
      Value sum = 0;
      for ( auto const& pos : pos_domain )
        {
          Point center;
          for ( size_t i = 0; i < dimension; ++i )
            center[i] = dx[i] * (pos[i] + 0.5);

          for ( auto const& point : domain )
            {
              Value value = 0.5 * ( 1. - std::tanh( 0.5*( (point-center).norm() - radius*dx[0] )/eps ));
              myImages[Label].setValue( point, value );
              sum += value;
            }

          ++Label;
        }
      
      return sum;
    }
  
  Value sumAllImages() const
    {
      Value sum = 0;
      //Domain const& domain = myImages[0].domain();
      for ( auto const& image : myImages )
        {
          for ( auto value : image )
            sum += value;
          //for ( auto const& point : domain )
          //  sum += image(point);
        }

      return sum;
    }

  Value sumOneImage( size_t aLabel ) const
    {
      Value sum = 0;
      Image const& image = myImages[aLabel];
      for ( auto value : image )
        sum += value;
      //Domain const& domain = image.domain();
      //for ( auto const& point : domain )
      //  sum += image(point);

      return sum;
    }

  Value sumSinAllImages(Value shift = 0) const
    {
      Value sum = 0;
      for ( auto const& image : myImages )
        {
          for ( auto value : image )
            sum += std::sin(value + shift);
        }

      return sum;
    }

  Value sumSinOneImage( size_t aLabel, Value shift = 0) const
    {
      Value sum = 0;
      Image const& image = myImages[aLabel];
      for ( auto value : image )
        sum += std::sin(value+shift);

      return sum;
    }

  size_t area() const
    {
      size_t total = 0;
      total += sizeof( std::vector<Image> );
      total += myImages.size() * sizeof(Image);
      total += myImages.size() * myImages[0].domain().size() * sizeof(Value);
      return total;
    }

private:
  std::vector<Image> myImages;
};

template <
  typename TMultiImage,
  size_t N
>
class BenchMultiImage
{
public:
  using MultiImage = TMultiImage;
  using Domain = typename MultiImage::Domain;
  using Point = typename Domain::Point;
  using Value = typename MultiImage::Value;
  using Approximation = typename MultiImage::Approximation;

  BenchMultiImage( Domain const& domain, Approximation const& approx )
    : myMultiImage{ domain, approx }
  {}

  inline size_t nImages() const
    {
      return ipow(N, Domain::dimension);
    }

  Value generateData( Value radius, Value eps )
    {
      const size_t dimension = Domain::dimension;
      Domain const& domain = myMultiImage.domain();

      // Shift between each phase center
      const Point extent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
      Value dx[dimension];
      for (size_t i = 0; i < dimension; ++i)
        dx[i] = Value(extent[i]) / N;

      const size_t L = ipow(N, dimension);
      Point center[L];
      Domain pos_domain{ Point::diagonal(0), Point::diagonal(N-1) };
      size_t Label = 0;
      for ( auto const& pos : pos_domain )
        {
          for (size_t i = 0; i < dimension; ++i)
              center[Label][i] = dx[i] * (pos[i] + 0.5);

          ++Label;
        }

      Value sum = 0;
      for ( auto const& point : domain )
        {
          for ( size_t i = 0; i < L; ++i )
            {
              Value value = 0.5 * ( 1. - std::tanh( 0.5*( (point-center[i]).norm() - radius*dx[0] )/eps ));
              myMultiImage.setValue( point, i, value );
              sum += value;
            }
        }

      return sum;
    }

  Value sumAllImages() const
    {
      Value sum = 0;
      Domain const& domain = myMultiImage.domain();
      for ( auto const& point : domain )
        for ( auto const& pair : myMultiImage(point) )
          sum += pair.second;

      return sum;
    }

  Value sumOneImage( size_t aLabel ) const
    {
      Value sum = 0;
      /*
      // Old fashion
      Domain const& domain = myMultiImage.getBoundingBox(aLabel);
      for ( auto const& point : domain )
        sum += myMultiImage.getValue(point, aLabel);
      */
      
      using ImageView = DGtal::ImageView<const MultiImage, DGtal::image_view::BoundingBoxAsDomain>;
      ImageView image_view( myMultiImage, aLabel );
      image_view.buffer() = Point::diagonal(0);

      for ( auto value : image_view )
        sum += value;

      return sum;
    }

  Value sumSinAllImages(Value shift = 0) const
    {
      Value sum = 0;
      Domain const& domain = myMultiImage.domain();
      for ( auto const& point : domain )
        for ( auto const& pair : myMultiImage(point) )
          sum += std::sin(pair.second+shift);

      return sum;
    }

  Value sumSinOneImage( size_t aLabel, Value shift = 0 ) const
    {
      Value sum = 0;
      /*
      // Old fashion
      Domain const& domain = myMultiImage.getBoundingBox(aLabel);
      for ( auto const& point : domain )
        sum += std::sin(myMultiImage.getValue(point, aLabel)+shift);
      */

      using ImageView = DGtal::ImageView<const MultiImage, DGtal::image_view::BoundingBoxAsDomain>;
      ImageView image_view( myMultiImage, aLabel );
      image_view.buffer() = Point::diagonal(0);

      for ( auto value : image_view )
        sum += std::sin(value + shift);

      return sum;
    }

  template <
    typename TData,
    unsigned int L,
    typename TWord,
    unsigned int NN,
    unsigned int M
  >
  size_t areaOfLabelledMap( LabelledMap<TData, L, TWord, NN, M> const& lmap ) const
    {
      size_t total = sizeof( LabelledMap<TData, L, TWord, NN, M> );
      const size_t size = lmap.size();
      if ( size > NN+1 )
        total += ( 1 + ( size - NN - 1 ) / M ) * ( M * sizeof( TData ) + sizeof(TData*) );
      return total;
    }

  size_t area () const
    {
      size_t total = 0;
      total += sizeof(myMultiImage);
      
      Domain const& domain = myMultiImage.domain();
      Point const extent = domain.upperBound() - domain.lowerBound() + Point::diagonal(1);
      size_t bb_size = sizeof(typename MultiImage::BoundingBox);
      for (size_t i = 0; i < Domain::dimension; ++i)
        {
          bb_size += extent[i] * sizeof(unsigned long);
        }
      total += ipow(N, Domain::dimension) * bb_size;

      size_t label_min = ipow(N, Domain::dimension);
      size_t label_max = 0;
      size_t label_sum = 0;
      for (auto const& point : domain)
        {
          auto && images = myMultiImage(point);
          total += areaOfLabelledMap( images );
          
          size_t nlabel = images.size();
          label_min = std::min( label_min, nlabel );
          label_max = std::max( label_max, nlabel );
          label_sum += nlabel;
        }

      std::cout << "\tNumber of label : min=" << label_min << " ; mean=" << double(label_sum)/domain.size() << " ; max=" << label_max << std::endl;

      return total;

    }

private:
  MultiImage myMultiImage;
};

template < typename TImages, typename T, typename ...Args >
void BenchIt ( std::string const& name, T radius, T eps, Args const& ... args )
{
  trace.beginBlock("----------- " + name + " ----------");
  trace.beginBlock("Allocating");
  TImages images( args... );
  trace.endBlock();
  
  trace.beginBlock("Initializing");
  std::cout << "\tsum = " << images.generateData(radius, eps) << std::endl;
  trace.endBlock();

  trace.beginBlock("Summing all images");
  std::cout << "\tsum = " << images.sumAllImages() << std::endl;
  trace.endBlock();

  trace.beginBlock("Summing all images one by one");
  T sum = 0;
  for (size_t i = 0; i < images.nImages(); ++i)
    sum += images.sumOneImage(i);
  std::cout << "\tsum = " << sum << std::endl;
  trace.endBlock();

  trace.beginBlock("Summing all images with sinus");
  std::cout << "\tsum = " << images.sumSinAllImages() << std::endl;
  trace.endBlock();

  trace.beginBlock("Summing all images with sinus, one by one");
  sum = 0;
  for (size_t i = 0; i < images.nImages(); ++i)
    sum += images.sumSinOneImage(i);
  std::cout << "\tsum = " << sum << std::endl;
  trace.endBlock();
  
  T const pi = boost::math::constants::pi<T>();
  trace.beginBlock("Summing all images with sinus+pi");
  std::cout << "\tsum = " << images.sumSinAllImages(pi) << std::endl;
  trace.endBlock();

  trace.beginBlock("Summing all images with sinus+pi, one by one");
  sum = 0;
  for (size_t i = 0; i < images.nImages(); ++i)
    sum += images.sumSinOneImage(i, pi);
  std::cout << "\tsum = " << sum << std::endl;
  trace.endBlock();

  trace.endBlock();

  std::cout << "\tMemory usage = " << double(images.area()/1024)/1024 << "Mo" << std::endl;
}


int main()
{
  using real = double;

  static const size_t D = 2;  ///< Space dimension.
  static const size_t N = 8;  ///< Number of images per dimension.
  static const size_t X = 1024; ///< Space size (in each dimension).
  static const size_t M = 5; ///< Additional capacity of LabelledMap

  const real radius = std::sqrt(2)/2; ///< Radius of the phase as ratio of the cell size.
  const real eps = 1; ///< Epsilon in phase-field initialization.

  static const size_t L = ipow(N, D); ///< Total number of images.

  using Space = SpaceND<D, DGtal::int32_t>;
  using Domain = HyperRectDomain<Space>;
  using Point = typename Domain::Point;

  using ImageContainerBySTLVector = ImageContainerBySTLVector< Domain, real >;


  using AABB = AxisAlignedBoundingBox< Domain, size_t>;
  using NoBB = NoBoundingBox<Domain>;

  using NoApprox  = approximations::NoValueApproximation<real>;
  using NegApprox = approximations::NegativeTolValueApproximation<real>;

  using LabelledMap1 = LabelledMap<real, L, size_t, 1, M>;
  using LabelledMap2 = LabelledMap<real, L, size_t, 2, M>;
  using LabelledMap3 = LabelledMap<real, L, size_t, 3, M>;
  using LabelledMap4 = LabelledMap<real, L, size_t, 4, M>;

  using ApproxMultiImage = ApproximatedMultiImage<Domain, LabelledMap4, NegApprox, NoBB>;

  Domain domain( Point::diagonal(0), Point::diagonal(X-1) );

  using std::cout; using std::endl;

  cout << "Benchmark in dimension " << D << " on a domain of size " << X << "^" << D << " with " << L << " images." << endl;
  cout << "Each image is initialized with the phase-field (eps=" << eps << ") corresponding to a ball of radius " << (X*radius/N) << endl;
  cout << endl;
  
  BenchIt< BenchVectorOfImages<ImageContainerBySTLVector, N> >("vector<ImageContainerBySTLVector>", radius, eps, domain); std::cout << std::endl;

  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap1, NoApprox, NoBB>, N> >("ApproximatedMultiImage - N=1 - no approx - no BB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap2, NoApprox, NoBB>, N> >("ApproximatedMultiImage - N=2 - no approx - no BB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap3, NoApprox, NoBB>, N> >("ApproximatedMultiImage - N=3 - no approx - no BB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap4, NoApprox, NoBB>, N> >("ApproximatedMultiImage - N=4 - no approx - no BB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;

  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap1, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=1 - approx 1e-10 - no BB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap2, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=2 - approx 1e-10 - no BB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap3, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=3 - approx 1e-10 - no BB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap4, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=4 - approx 1e-10 - no BB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap1, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=1 - approx 1e-4 - no BB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap2, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=2 - approx 1e-4 - no BB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap3, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=3 - approx 1e-4 - no BB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap4, NegApprox, NoBB>, N> >("ApproximatedMultiImage - N=4 - approx 1e-4 - no BB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap1, NoApprox, AABB>, N> >("ApproximatedMultiImage - N=1 - no approx - AABB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap2, NoApprox, AABB>, N> >("ApproximatedMultiImage - N=2 - no approx - AABB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap3, NoApprox, AABB>, N> >("ApproximatedMultiImage - N=3 - no approx - AABB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap4, NoApprox, AABB>, N> >("ApproximatedMultiImage - N=4 - no approx - AABB", radius, eps, domain, NoApprox{} ); std::cout << std::endl;
  
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap1, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=1 - approx 1e-10 - AABB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap2, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=2 - approx 1e-10 - AABB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap3, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=3 - approx 1e-10 - AABB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap4, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=4 - approx 1e-10 - AABB", radius, eps, domain, NegApprox{1e-10} ); std::cout << std::endl;

  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap1, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=1 - approx 1e-4 - AABB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap2, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=2 - approx 1e-4 - AABB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap3, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=3 - approx 1e-4 - AABB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;
  BenchIt< BenchMultiImage< ApproximatedMultiImage<Domain, LabelledMap4, NegApprox, AABB>, N> >("ApproximatedMultiImage - N=4 - approx 1e-4 - AABB", radius, eps, domain, NegApprox{1e-4} ); std::cout << std::endl;

  BenchIt< BenchVectorOfImages<ImageContainerBySTLVector, N> >("vector<ImageContainerBySTLVector>", radius, eps, domain); std::cout << std::endl;

  return 0;
}
/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

