#include <mplapack_mpfr.h>

int main() {
    mpfr_class value(1);
    return Risnan(value) ? 1 : 0;
}
