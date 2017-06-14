#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "ImageTransform.h"
#include "ComponentWiseImageOperator.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace DGtal;
using namespace Z2i;

void ImageFill( ImageContainerBySTLVector<Domain, double> & anImage, double aValue )
{
  ImageTransform( anImage, anImage, [aValue] ( Point const&, double ) { return aValue; } );
}

void ImageDouble(
    ImageContainerBySTLVector<Domain, double> & anInputImage,
    ImageContainerBySTLVector<Domain, double> & anOutputImage )
{
  using Image = ImageContainerBySTLVector<Domain, double>;
  using Context = ComponentWiseImageContext<Image>;
  auto const imageOp = ComponentWiseImageOperator(
          [] ( Context::Value const& value )
          {
            return 2*value();
          }
  );

  imageOp( anInputImage ).saveTo( anOutputImage );
}

int main( int argc, char* argv[] )
{

  using real = double;
  using Image = ImageContainerBySTLVector<Domain, real>;
  const double a = std::stod( argv[1] );

  Domain domain( {0,0}, {10,10} );
  Image image1( domain );
  Image image2( domain );

  const Point pt( 3, 7 );

  std::cout << image1( pt ) << std::endl;
  ImageFill( image1, a );
  std::cout << image1( pt ) << std::endl;

  ImageTransform( image1, image2, [] ( Point const& p, real v ) { return p[0] == 3 ? v+1 : v-1; } );
  std::cout << image2( pt ) << std::endl;
  std::cout << image2( {4, 7} ) << std::endl;

  ImageDouble( image1, image2 );
  std::cout << image2( pt ) << std::endl;

  using Context = ComponentWiseImageContext<Image>;

  image1
      >> ComponentWiseImageOperator(
              [a] ( Context::Value const& value )
              {
                return value() * a;
              })
      >> image2;
  std::cout << image2( pt ) << std::endl;

  return 0;
}
