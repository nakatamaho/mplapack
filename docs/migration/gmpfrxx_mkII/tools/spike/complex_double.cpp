#include <gmpfrxx_mkII.h>

#include <complex>

int main() {
    std::complex<double> source(1.25, -2.5);
    mpfrxx::mpc_class target(source);
    return mpfrxx::real(target) == 1.25 && mpfrxx::imag(target) == -2.5 ? 0
                                                                        : 1;
}
