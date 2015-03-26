#include <iostream>

#include "ValueApproximations.h"

int main()
{
    using namespace std;

    using real = double;
    
    BOOST_CONCEPT_ASSERT( (DGtal::concepts::CValueApproximation< DGtal::approximation::NoValueApproximation<real> >) );
    BOOST_CONCEPT_ASSERT( (DGtal::concepts::CValueApproximation< DGtal::approximation::ZeroValueApproximation<real> >) );
    BOOST_CONCEPT_ASSERT( (DGtal::concepts::CValueApproximation< DGtal::approximation::ZeroTolValueApproximation<real> >) );
    BOOST_CONCEPT_ASSERT( (DGtal::concepts::CValueApproximation< DGtal::approximation::NegativeValueApproximation<real> >) );
    BOOST_CONCEPT_ASSERT( (DGtal::concepts::CValueApproximation< DGtal::approximation::NegativeTolValueApproximation<real> >) );
    
    {
        cout << "NoValueApproximation:" << endl;
        const auto approx = DGtal::approximations::NoValueApproximation<real>{};
        cout << ( approx.eval(1) == false ? "OK" : "KO" ) << endl;
    }

    {
        cout << "ZeroValueApproximation:" << endl;
        const auto approx = DGtal::approximations::ZeroValueApproximation<real>{};
        cout << ( !approx.eval(-1) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(0) ? "OK" : "KO" ) << endl;
    }

    {
        cout << "ZeroTolValueApproximation:" << endl;
        const auto approx = DGtal::approximations::ZeroTolValueApproximation<real>{1e-10};
        cout << ( approx.eval(0) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(1e-15) ? "OK" : "KO" ) << endl;
        cout << ( !approx.eval(1) ? "OK" : "KO" ) << endl;
    }

    {
        cout << "NegativeTolValueApproximation:" << endl;
        const auto approx = DGtal::approximations::NegativeTolValueApproximation<real>{1e-10};
        cout << ( approx.eval(-1) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(0) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(1e-15) ? "OK" : "KO" ) << endl;
        cout << ( !approx.eval(1) ? "OK" : "KO" ) << endl;
    }


    return 0;
}
