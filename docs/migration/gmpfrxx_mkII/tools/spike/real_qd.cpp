#include <gmpfrxx_mkII.h>
#include "mpblas_qd.h"

int main() {
    qd_real source = qd_real(1.0) + ldexp(qd_real(1.0), -100);
    mpfrxx::mpfr_class target(source);
    return target - mpfrxx::mpfr_class(1.0) > 0.0 ? 0 : 1;
}
