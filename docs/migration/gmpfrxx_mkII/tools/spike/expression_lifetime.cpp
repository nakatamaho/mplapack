#include <gmpfrxx_mkII.h>

int main() {
    mpfrxx::mpfr_class result =
        (mpfrxx::mpfr_class(1.0) + mpfrxx::mpfr_class(2.0)) *
        (mpfrxx::mpfr_class(3.0) + mpfrxx::mpfr_class(4.0));
    mpfrxx::mpc_class complex_result =
        (mpfrxx::mpc_class(1.0, 2.0) + mpfrxx::mpc_class(3.0, 4.0)) *
        mpfrxx::mpc_class(2.0, 0.0);
    return result == 21.0 && mpfrxx::real(complex_result) == 8.0 &&
                   mpfrxx::imag(complex_result) == 12.0
               ? 0
               : 1;
}
