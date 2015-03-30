#include <iostream>
#include <cstddef>

#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

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
    using Approximation = DGtal::approximations::ZeroTolValueApproximation<real>;
    using Space = SpaceND<N, int>;
    using Domain = HyperRectDomain<Space>;
    using Point = typename Domain::Point;
    using BoundingBox = AxisAlignedBoundingBox< Domain, unsigned int>;

    //using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation>;
    using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation, BoundingBox>;

    const Domain domain({0,0}, {1000,100,100});
    const Approximation approx{1};
    ApproximatedMultiImage images(domain, approx);

    cout << images.getValue( {100,100,100}, 0 ) << endl;
    images.setValue( {100,100,100}, 0, 1.1 );
    cout << images.getValue( {100,100,100}, 0 ) << endl;

    using ImageView = ImageView<ApproximatedMultiImage, image_view::BoundingBoxAsDomain>;
    //using ImageView = ImageView<ApproximatedMultiImage>;
    ImageView image_view{ images, 0 };
    image_view.buffer() = Domain::Point::diagonal(1);

    cout << image_view.domain() << endl;

    cout << images[0]({100,100,100}) << endl;

    auto const& cimages = images;
    cout << cimages[0].domain() << endl;
    cout << cimages[0]({100,100,100}) << endl;
    cout << cimages.getBBImage(0).domain() << endl;


    return 0;
}
