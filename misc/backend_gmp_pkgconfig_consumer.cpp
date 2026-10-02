#include <mplapack_gmp.h>

int main() {
    mpf_class value(1);
    return Risnan(value) ? 1 : 0;
}
