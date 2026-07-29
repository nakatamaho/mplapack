#include <gmpfrxx_mkII.h>
#include "mplapack_config.h"
#include "mplapack_utils_gmp.h"

int main() {
    gmpxx::mpf_class value = 2.0;
    char buffer[1024];
    auto square = pow2(value);
    auto signed_value = sign(value, -value);
    auto nearest = nint(value);
    auto circle = pi(value);
    sprintnum(buffer, value);
    return square > value && signed_value < 0 && nearest == 2 && circle > value;
}
