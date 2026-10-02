#include <mplapack_mpfr.h>

int main() {
    mpfr::mpreal value(1);
    return Risnan(value) ? 1 : 0;
}
