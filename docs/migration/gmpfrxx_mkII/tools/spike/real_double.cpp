#include <gmpfrxx_mkII.h>

int main() {
    double source = 1.25;
    mpfrxx::mpfr_class target(source);
    return target == 1.25 ? 0 : 1;
}
