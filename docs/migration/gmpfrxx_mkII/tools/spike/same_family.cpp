#include <gmpfrxx_mkII.h>

int main() {
    mpfrxx::mpfr_class a = 2.0;
    mpfrxx::mpfr_class b = 3.0;
    mpfrxx::mpfr_class real_result = a * b + sqrt(a);
    mpfrxx::mpc_class z(a, b);
    mpfrxx::mpc_class complex_result = z * z + z;
    gmpxx::mpf_class f = 2.0;
    gmpxx::mpf_class g = 3.0;
    gmpxx::mpf_class gmp_result = f * g + f;
    gmpxx::mpfc_class c(f, g);
    gmpxx::mpfc_class gmp_complex_result = c * c + c;
    return real_result > a && mpfrxx::real(complex_result) != a &&
                   gmp_result > f && gmpxx::real(gmp_complex_result) != f
               ? 0
               : 1;
}
