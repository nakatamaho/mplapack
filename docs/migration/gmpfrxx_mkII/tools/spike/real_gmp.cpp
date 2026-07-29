#include <gmpfrxx_mkII.h>

int main() {
    gmpxx::mpf_class source("1.25");
    mpfrxx::mpfr_class target(source);
    return target == 1.25 ? 0 : 1;
}
