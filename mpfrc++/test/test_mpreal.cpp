#include "mpreal.h"
#include <iostream>

using namespace mpfr;
using namespace std;

int main(int argc, char *argv[]) {

    printf("%ld\n", mpreal::get_default_prec());
    mpreal::set_default_prec(128);
    printf("%ld\n", mpreal::get_default_prec());

    return 0;
}
