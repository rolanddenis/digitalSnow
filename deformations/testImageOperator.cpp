#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "ComponentWiseImageOperator.h"

#include <iostream>
#include <string>
#include <cmath>
#include <boost/core/demangle.hpp>
#include <typeinfo>

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

struct Functor
{
  static std::size_t cnt;
  Functor() { ++cnt; }
  Functor( Functor const& ) { ++cnt; }
  Functor( Functor && )     { ++cnt; }
  ~Functor() { --cnt; }

  template < typename TContext >
  double operator() ( TContext const& ct ) const
    {
      return std::cos(ct.value);
    }
};

std::size_t Functor::cnt = 0;

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
  std::cout << "type(imageOp) = " << boost::core::demangle( typeid(imageOp).name() ) << std::endl;
  auto const result = imageOp( image );
  std::cout << "type(result) = " << boost::core::demangle( typeid(result).name() ) << std::endl;

  auto const result_rvalue = makeComponentWiseImageOperator( 
          [] ( Context ct )
          {
            return std::cos(ct.value);
          }
  )( image );
  std::cout << "type(result_rvalue) = " << boost::core::demangle( typeid(result_rvalue).name() ) << std::endl;

  auto const fn = dummy2<Context>;
  auto const result2 = makeComponentWiseImageOperator( fn )( image );
  std::cout << "type(result2) = " << boost::core::demangle( typeid(result2).name() ) << std::endl;

  auto const result3 = makeComponentWiseImageOperator( dummy2<Context> )( image );
  std::cout << "type(result3) = " << boost::core::demangle( typeid(result3).name() ) << std::endl;

  std::cout << "cnt = " << Functor::cnt << std::endl;
  auto const result4 = makeComponentWiseImageOperator( Functor{} )( image );
  std::cout << "type(result4) = " << boost::core::demangle( typeid(result4).name() ) << std::endl;
  std::cout << "cnt = " << Functor::cnt << std::endl;

  auto const fn2 = Functor{};
  std::cout << "cnt = " << Functor::cnt << std::endl;
  auto const imageOp2 = makeComponentWiseImageOperator( fn2 );
  std::cout << "cnt = " << Functor::cnt << std::endl;
  auto const result5 = imageOp2( image );
  std::cout << "cnt = " << Functor::cnt << std::endl;
  std::cout << "type(imageOp2) = " << boost::core::demangle( typeid(imageOp2).name() ) << std::endl;
  std::cout << "type(result5) = " << boost::core::demangle( typeid(result5).name() ) << std::endl;

  const Point pt(0, 0);

  image.setValue(pt, std::stod( argv[1] ) );
  std::cout << image(pt) << std::endl;
  std::cout << result(pt) << std::endl;
  std::cout << result_rvalue(pt) << std::endl;
  std::cout << dummy( image, pt ) << std::endl;
  std::cout << result2(pt) << std::endl;
  std::cout << result3(pt) << std::endl;
  std::cout << result4(pt) << std::endl;
  return 0;
}
