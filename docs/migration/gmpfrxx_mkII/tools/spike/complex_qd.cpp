#include <gmpfrxx_mkII.h>
#include "mpblas_qd.h"

int main() {
    qd_real real = qd_real(1.0) + ldexp(qd_real(1.0), -100);
    qd_complex source(real, qd_real(-2.0));
    mpfrxx::mpc_class target(source);
    return mpfrxx::real(target) - mpfrxx::mpfr_class(1.0) > 0.0 &&
                   mpfrxx::imag(target) == -2.0
               ? 0
               : 1;
}
