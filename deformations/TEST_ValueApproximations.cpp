#include <iostream>

#include "ValueApproximations.h"

int main()
{
    using namespace std;

    using real = double;
    
    {
        cout << "NoValueApproximation:" << endl;
        const auto approx = DGtal::approximation::NoValueApproximation<real>{};
        cout << ( approx.eval(1) == false ? "OK" : "KO" ) << endl;
    }

    {
        cout << "ZeroValueApproximation:" << endl;
        const auto approx = DGtal::approximation::ZeroValueApproximation<real>{};
        cout << ( !approx.eval(-1) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(0) ? "OK" : "KO" ) << endl;
    }

    {
        cout << "ZeroTolValueApproximation:" << endl;
        const auto approx = DGtal::approximation::ZeroTolValueApproximation<real>{1e-10};
        cout << ( approx.eval(0) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(1e-15) ? "OK" : "KO" ) << endl;
        cout << ( !approx.eval(1) ? "OK" : "KO" ) << endl;
    }

    {
        cout << "NegativeTolValueApproximation:" << endl;
        const auto approx = DGtal::approximation::NegativeTolValueApproximation<real>{1e-10};
        cout << ( approx.eval(-1) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(0) ? "OK" : "KO" ) << endl;
        cout << ( approx.eval(1e-15) ? "OK" : "KO" ) << endl;
        cout << ( !approx.eval(1) ? "OK" : "KO" ) << endl;
    }


    return 0;
}
