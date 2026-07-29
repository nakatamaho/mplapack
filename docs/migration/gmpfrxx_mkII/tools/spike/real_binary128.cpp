#include <gmpfrxx_mkII.h>
#include "mpblas_binary128.h"

int main() {
    mplapack_binary128_t source =
        static_cast<mplapack_binary128_t>(1.0) + 0x1p-100F128;
    mpfrxx::mpfr_class target(source);
    return target - mpfrxx::mpfr_class(1.0) > 0.0 ? 0 : 1;
}
