#include <stdio.h>
#include <gmpxx.h>
#include <mpc_class.h>
#include <complex>

int main(int argc, char *argv[])
{
  mpc_class a;
  std::complex<double> b;
  a.real() = 100.0;
  a.imag() = 3.0;
  b = 1.0;
  b.imag() = 3.0;
  gmp_printf("%10.15Fe + %10.15Fe i \n", a.real().get_mpf_t(), a.imag().get_mpf_t());
  printf("%10.15e + %10.15e i \n", b.real(), b.imag());
  return 1;
}
