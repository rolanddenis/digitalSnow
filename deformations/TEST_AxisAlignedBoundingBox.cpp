#include <iostream>

#include <DGtal/base/Common.h>
#include <DGtal/kernel/SpaceND.h>
#include "AxisAlignedBoundingBox.h"

int main()
{
    using namespace DGtal;

    using Space = SpaceND<2, int>;
    using Domain = HyperRectDomain<Space>;
    using BoundingBox = AxisAlignedBoundingBox< Domain, unsigned long, 0>;
    using DefaultBoundingBox = AxisAlignedBoundingBox< Domain, unsigned long>;

    Domain domain({1,2}, {10,20});
    BoundingBox toto(domain);

    using namespace std;

    cout << sizeof(BoundingBox) << " =? " << sizeof(DefaultBoundingBox) << endl;
    cout << sizeof(std::vector<unsigned long>) << endl;
    cout << sizeof(Domain) << endl;
    cout << sizeof(typename Space::Point) << endl;

    cout << toto << endl;
    toto.addPoint( { 5, 6 } );
    cout << toto << endl;
    toto.addPoint( { 7, 6 } );
    cout << toto << endl;
    toto.addPoint( { 7, 3 } );
    cout << toto << endl;
    toto.removePoint( { 7, 6 } );
    cout << toto << endl;
    toto.removePoint( { 5, 6 } );
    cout << toto << endl;
    cout << toto.getBoundingBox( Space::Point::diagonal(2) ) << endl;
    toto.removePoint( { 7, 3 } );
    cout << toto << endl;

    return 0;
}
