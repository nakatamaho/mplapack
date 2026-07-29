#include <gmpfrxx_mkII.h>
#include "mpblas_dd.h"

int main() {
    dd_real real = dd_real(1.0) + ldexp(dd_real(1.0), -100);
    dd_complex source(real, dd_real(-2.0));
    mpfrxx::mpc_class target(source);
    return mpfrxx::real(target) - mpfrxx::mpfr_class(1.0) > 0.0 &&
                   mpfrxx::imag(target) == -2.0
               ? 0
               : 1;
}
