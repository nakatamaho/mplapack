#define ___MPLAPACK_MPLAPACK_INIT___
#include "mpreal.h"
#include <iostream>

using namespace mpfr;
using namespace std;

int main(int argc, char *argv[]) {
    mp_rnd_t RND = mpreal::default_rnd;
    mpreal::default_rnd = RND;
    mpreal::default_prec = 256;
    mpreal source = "0.1";
    source.set_prec(128);
    mpreal::default_prec = 64;
    mpreal copied(source);
    if (copied.get_prec() != source.get_prec()) {
        cout << "copy constructor precision test failed" << endl;
        return 1;
    }
    if (copied != source) {
        cout << "copy constructor value test failed" << endl;
        return 1;
    }

    mpreal::default_prec = 256;
    mpreal a = "0.1";
    a.set_prec(8);

    mpreal b = "0.1";
    a.set_prec(128);
    mpreal::default_prec = 256;
    mpreal d = "0.2";

    mpreal c = a + b;
    mpfr_printf("%50.100Re\n", a);
    mpfr_printf("%50.100Re\n", b);
    mpfr_printf("%50.100Re\n", c);
    mpfr_printf("%50.100Re\n", d);

    mpfr_printf("  123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789\n");    
    mpfr_printf(" %50.130Rb\n", a);
    mpfr_printf(" %50.130Rb\n", b);
    mpfr_printf("%50.130Rb\n", c);
    mpfr_printf("%50.130Rb\n", d);

    return 0;
}
