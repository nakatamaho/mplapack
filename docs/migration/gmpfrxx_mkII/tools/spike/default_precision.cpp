#include <gmpfrxx_mkII.h>

int main() {
    mpfrxx::mpfr_class mpfr_value;
    gmpxx::mpf_class mpf_value;
    return mpfr_value.get_prec() == 512 && mpf_value.get_prec() == 512 ? 0 : 1;
}
