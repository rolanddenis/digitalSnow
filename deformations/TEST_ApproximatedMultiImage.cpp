#include <iostream>
#include <cstddef>

#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/CImage.h>
#include <DGtal/images/CConstImage.h>
#include <DGtal/images/CImage.h>

#include <boost/range/adaptor/reversed.hpp>

#include "AxisAlignedBoundingBox.h"
#include "ValueApproximations.h"
#include "ApproximatedMultiImage.h"
#include "ImageView.h"

int main()
{
    const std::size_t N = 3;

    using namespace DGtal;
    using namespace std;

    using real = double;
    using LabelledMap = DGtal::LabelledMap<real, 32, unsigned int, 2, 2>;
    //using Approximation = DGtal::approximations::ZeroValueApproximation<real>;
    //using Approximation = DGtal::approximations::ZeroTolValueApproximation<real>;
    using Approximation = DGtal::approximations::NegativeTolValueApproximation<real>;

    using Space = SpaceND<N, int>;
    using Domain = HyperRectDomain<Space>;
    using Point = typename Domain::Point;
    using BoundingBox = AxisAlignedBoundingBox< Domain, unsigned int>;

    //using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation>;
    using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation, BoundingBox>;
    using ImageView = ImageView<ApproximatedMultiImage, image_view::BoundingBoxAsDomain>;

    const Domain domain({0,0}, {1000,100,100});
    const Approximation approx{1};
    ApproximatedMultiImage images(domain, approx);
    

    cout << images.getValue( {100,100,100}, 0 ) << endl;
    images.setValue( {100,100,100}, 0, 1.1 );
    cout << images.getValue( {100,100,100}, 0 ) << endl;

    //using ImageView = ImageView<ApproximatedMultiImage>;
    ImageView image_view{ images, 0 };
    image_view.buffer() = Domain::Point::diagonal(1);

    cout << image_view.domain() << endl;

    real sum = 0;
    /*
    for ( auto value : image_view )
        sum += value;
    cout << "sum = " << sum << endl;
    */
    
    sum = 0;
    //for ( auto value : boost::adaptors::reverse(image_view) )
    //    sum += value;
    cout << image_view.domain() << endl;
    auto pt = image_view.begin();
    *pt = 2.;
    cout << image_view.domain() << endl;
    cout << *pt << endl;

    for ( auto it = image_view.rbegin(), it_end = image_view.rend(); it != it_end; ++it )
        sum += *it;
    cout << "sum = " << sum << endl;
    
    sum = 0;
    ImageView image_view1{ images, 1}; 
    image_view1.buffer() = Domain::Point::diagonal(0);
    for ( auto value : image_view1 )
        sum += value;
    cout << "sum = " << sum << endl;
    cout << image_view1.domain() << " : Iterators ? " << ( image_view1.begin() == image_view1.end() ? "OK" : "KO" ) << endl;

    cout << images[0]({100,100,100}) << endl;

    BOOST_CONCEPT_ASSERT( (DGtal::concepts::CImage<ImageView>) );

    auto const& cimages = images;
    cout << cimages[0].domain() << endl;
    cout << cimages[0]({100,100,100}) << endl;
    cout << cimages.getBBImage(0).domain() << endl;

    ImageContainerBySTLVector<Domain, real> plain_image ( images.getBBImage(0) );
    cout << plain_image({100,100,100}) << endl;

    return 0;
}
