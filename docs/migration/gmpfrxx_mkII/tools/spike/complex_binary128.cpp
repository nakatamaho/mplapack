#include <gmpfrxx_mkII.h>
#include "mpblas_binary128.h"

#include <complex>

int main() {
    using source_real = mplapack_binary128_t;
    std::complex<source_real> source(
        static_cast<source_real>(1.0) + 0x1p-100F128,
        static_cast<source_real>(-2.0));
    mpfrxx::mpc_class target(source);
    return mpfrxx::real(target) - mpfrxx::mpfr_class(1.0) > 0.0 &&
                   mpfrxx::imag(target) == -2.0
               ? 0
               : 1;
}
