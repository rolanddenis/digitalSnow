#include <iostream>
#include <cstddef>

#include <DGtal/base/LabelledMap.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

#include "AxisAlignedBoundingBox.h"
#include "ValueApproximations.h"
#include "ApproximatedMultiImage.h"

int main()
{
    const std::size_t N = 2;

    using namespace DGtal;
    using namespace std;

    using real = double;
    using LabelledMap = DGtal::LabelledMap<real, 32, unsigned int, 2, 2>;
    //using Approximation = DGtal::approximations::ZeroValueApproximation<real>;
    using Approximation = DGtal::approximations::ZeroTolValueApproximation<real>;
    using Space = SpaceND<N, int>;
    using Domain = HyperRectDomain<Space>;

    using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation>;

    const Domain domain({1,2}, {10,20});
    const Approximation approx{1};
    ApproximatedMultiImage images(domain, approx);

    cout << images.getValue( {2,2}, 0 ) << endl;
    images.setValue( {2,2}, 0, 1.1 );
    cout << images.getValue( {2,2}, 0 ) << endl;


    return 0;
}
