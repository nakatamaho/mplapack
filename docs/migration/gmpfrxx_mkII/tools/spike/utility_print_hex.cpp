#include <gmpfrxx_mkII.h>

#include <iomanip>
#include <sstream>
#include <string>

int main() {
    mpfrxx::mpfr_class value = 2.0;
    mpfrxx::mpfr_class root = sqrt(value);
    mpfrxx::mpc_class complex_value(value, root);
    mpfrxx::mpfr_class magnitude = abs(complex_value);
    std::ostringstream decimal;
    std::ostringstream hexadecimal;
    decimal << std::setprecision(40) << magnitude;
    hexadecimal << std::hexfloat << magnitude;
    return decimal.str().empty() || hexadecimal.str().empty() ? 1 : 0;
}
