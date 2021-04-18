#include "mpreal.h"
#include <iostream>

using namespace mpfr;
using namespace std;

// MPFR  C++ - version
mpreal mpfr_schwefel(mpreal x) { return "418.9829" - x * sin(sqrt(fabs(x))); }

int main(int argc, char *argv[]) {
    mpreal::set_default_prec(128);
    mpreal x = "420.968746359948568";

    cout.precision(15);
    cout << "Global minimum of the Schwefel function is = " << mpfr_schwefel(x) << endl;

    cout << "\nmachine epsilon = " << machine_epsilon(x.get_prec());
    cout << "\nmin = " << mpreal_min(x.get_prec());
    cout << "\nmax = " << mpreal_max(x.get_prec()) << endl;

    return 0;
}
