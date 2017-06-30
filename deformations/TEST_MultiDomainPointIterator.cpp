#include <iostream>
#include <cstddef>
#include <vector>

#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/kernel/domains/CDomain.h>

#include "MultiDomainPointIterator.h"

int main()
{
    const std::size_t N = 2;

    using namespace DGtal;
    using namespace std;

    using real      = double;
    using Space     = SpaceND<N, int>;
    using Domain    = HyperRectDomain<Space>;
    using Point     = typename Domain::Point;

    std::vector< Domain > multiDomain = {
        { Point{0,0}, Point{2,2} },
        { Point{3,2}, Point{6,2} },
        { Point{0,3}, Point{3,4} }
    };

    const auto domain_begin  = multiDomain.begin(); 
    const auto domain_end    = multiDomain.end(); 

    auto mdomain_it = MultiDomainPointIterator<decltype(domain_begin)>( domain_begin, domain_end );
    auto mdomain_end_it = MultiDomainPointIterator<decltype(domain_begin)>( domain_begin, domain_end, true );

    for ( ; mdomain_it != mdomain_end_it; ++mdomain_it )
        std::cout << *mdomain_it << std::endl;

    return 0;
}
