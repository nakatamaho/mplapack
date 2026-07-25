#include <gmpfrxx_mkII.h>
#include "mpc_class.h"

int main() {
    ::mpc_class legacy(1.0, 2.0);
    mpfrxx::mpc_class replacement(1.0, 2.0);
    gmpxx::mpfc_class replacement_gmp(1.0, 2.0);
    return legacy.real() == 1.0 && mpfrxx::real(replacement) == 1.0 &&
                   gmpxx::real(replacement_gmp) == 1.0
               ? 0
               : 1;
}
