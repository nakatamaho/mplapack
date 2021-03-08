#include "mpreal.h"
#include <iostream>

using namespace mpfr;
using namespace std;

int main(int argc, char *argv[]) {
    mp_rnd_t RND = mpreal::get_default_rnd();
    mpreal::set_default_rnd(RND);
    mpreal::set_default_prec(256);
    mpreal a = "0.1";
    a.set_prec(8);

    mpreal b = "0.1";
    a.set_prec(128);
    mpreal::set_default_prec(256);
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
