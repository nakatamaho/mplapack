#include <gmpfrxx_mkII.h>
#include "mpblas_binary80.h"

int main() {
    mplapack_binary80_t source =
        static_cast<mplapack_binary80_t>(1.0) + 0x1p-60F64x;
    mpfrxx::mpfr_class target(source);
    return target - mpfrxx::mpfr_class(1.0) > 0.0 ? 0 : 1;
}
