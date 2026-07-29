#include <gmpfrxx_mkII.h>
#include "mpblas_dd.h"

int main() {
    dd_real source = dd_real(1.0) + ldexp(dd_real(1.0), -100);
    mpfrxx::mpfr_class target(source);
    return target - mpfrxx::mpfr_class(1.0) > 0.0 ? 0 : 1;
}
