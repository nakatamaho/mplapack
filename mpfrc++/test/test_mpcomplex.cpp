#include "mpreal.h"
#include "mpcomplex.h"
#include <iostream>

using namespace mpfr;
using namespace std;

int main(int argc, char *argv[]) {
    printf("%ld\n", mpcomplex::get_default_real_prec());
    return 0;
}
