#include "mpreal.h"
#include <iostream>

using namespace mpfr;
using namespace std;

#define buflen 1024
int main(int argc, char *argv[]) {
    mpreal::set_default_prec(256);
    mp_rnd_t RND = mpreal::get_default_rnd();
    mpreal::set_default_rnd(RND);
    const char *PI = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095";
    const char *PI_256 = "1.1001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111000001110011010001001p+1";
    const char *PI_128 = "1.1001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111000001110011010001000p+1";
    const char *PI_64 = "1.1001001000011111101101010100010001000010110100011000010001101010000000000000000000000000000000000000000000000000000000000000000000p+1";
    const char *PI_double = "1.1001001000011111101101010100010001000010110100011000000000000000000000000000000000000000000000000000000000000000000000000000000000p+1";
    const char *PI_float = "1.1001001000011111101101100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000p+1";
    char PI_256_TEST[buflen];
    char PI_128_TEST[buflen];
    char PI_64_TEST[buflen];
    char PI_double_TEST[buflen];
    char PI_float_TEST[buflen];

    mpreal pi = PI;
    mpreal pi_tranc128 = pi;
    pi_tranc128.set_prec(128);
    mpreal pi_tranc64 = pi;
    pi_tranc64.set_prec(64);
    double pi_double = (double)pi;
    double pi_float = (float)pi;

    printf("pi and trancation\n");
    mpfr_printf("%s\n", PI);
    mpfr_printf("%50.100Re\n", pi);
    mpfr_printf("%50.100Re\n", pi_tranc128);
    mpfr_printf("%50.100Re\n", pi_tranc64);
    mpfr_printf("%50.100Re\n", mpreal(pi_double));
    mpfr_printf("%50.100Re\n", mpreal(pi_float));

    mpfr_printf("  123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789\n");
    mpfr_printf("%50.130Rb\n", pi);
    mpfr_printf("%50.130Rb\n", pi_tranc128);
    mpfr_printf("%50.130Rb\n", pi_tranc64);
    mpfr_printf("%50.130Rb\n", mpreal(pi_double));
    mpfr_printf("%50.130Rb\n", mpreal(pi_float));

    mpfr_snprintf(PI_256_TEST, buflen, "%50.130Rb", pi);
    mpfr_snprintf(PI_128_TEST, buflen, "%50.130Rb", pi_tranc128);
    mpfr_snprintf(PI_64_TEST, buflen, "%50.130Rb", pi_tranc64);
    mpfr_snprintf(PI_double_TEST, buflen, "%50.130Rb", mpreal(pi_double));
    mpfr_snprintf(PI_float_TEST, buflen, "%50.130Rb", mpreal(pi_float));

    if (strncmp(PI_256_TEST, PI_256, buflen) != 0) {
        printf("ng\n");
        printf("%s\n", PI_256_TEST);
        printf("%s\n", PI_256);
    } else {
        printf("ok\n");
    }

    if (strncmp(PI_128_TEST, PI_128, buflen) != 0) {
        printf("ng\n");
        printf("%s\n", PI_128_TEST);
        printf("%s\n", PI_128);
    } else {
        printf("ok\n");
    }

    if (strncmp(PI_64_TEST, PI_64, buflen) != 0) {
        printf("ng\n");
        printf("%s\n", PI_64_TEST);
        printf("%s\n", PI_64);
    } else {
        printf("ok\n");
    }

    if (strncmp(PI_double_TEST, PI_double, buflen) != 0) {
        printf("ng\n");
        printf("%s\n", PI_double_TEST);
        printf("%s\n", PI_double);
    } else {
        printf("ok\n");
    }

    if (strncmp(PI_float_TEST, PI_float, buflen) != 0) {
        printf("ng\n");
        printf("%s\n", PI_float_TEST);
        printf("%s\n", PI_float);
    } else {
        printf("ok\n");
    }
    return 0;
}
