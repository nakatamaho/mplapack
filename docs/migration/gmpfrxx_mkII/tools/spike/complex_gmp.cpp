#include <gmpfrxx_mkII.h>

int main() {
    gmpxx::mpfc_class source(gmpxx::mpf_class("1.25"),
                             gmpxx::mpf_class("-2.5"));
    mpfrxx::mpc_class target(source);
    return mpfrxx::real(target) == 1.25 && mpfrxx::imag(target) == -2.5 ? 0
                                                                        : 1;
}
