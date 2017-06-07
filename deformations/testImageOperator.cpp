#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "ComponentWiseImageOperator.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace DGtal;
using namespace Z2i;

template < typename TImage, typename TPoint >
double dummy( TImage const& anImage, TPoint const& aPoint )
{
  using Context = ComponentWiseImageContext<typename TImage::Point, typename TImage::Value>;
  return makeComponentWiseImageOperator( 
          [] ( Context ct )
          {
            return std::cos(ct.value);
          }
  )( anImage )( aPoint );
}

template < typename TContext >
double dummy2( TContext const& context )
{
  return std::cos( context.value );
}

int main( int argc, char* argv[] )
{
  using real = double;
  using Image = ImageContainerBySTLVector<Domain, real>;
  using Context = ComponentWiseImageContext<Image::Point, real>;

  Domain domain( {0,0}, {10,10} );
  Image image( domain );

  auto const imageOp = makeComponentWiseImageOperator( 
          [] ( Context ct )
          {
            return 2*ct.value;
          }
  );
  auto const result = imageOp( image );

  auto const result_rvalue = makeComponentWiseImageOperator( 
          [] ( Context ct )
          {
            return std::cos(ct.value);
          }
  )( image );

  auto const fn = dummy2<Context>;
  auto const result2 = makeComponentWiseImageOperator( fn )( image );

  const Point pt(0, 0);

  image.setValue(pt, std::stod( argv[1] ) );
  std::cout << image(pt) << std::endl;
  std::cout << result(pt) << std::endl;
  std::cout << result_rvalue(pt) << std::endl;
  std::cout << dummy( image, pt ) << std::endl;
  return 0;
}
