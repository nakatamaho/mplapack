#include <stdio.h>
#include <mpc.h>

int main()
{
  mpfr_t x;
  mpc_t y;

  mpfr_init2(x, 100);
  mpc_init2(y, 100);

  mpc_clear(y);
  mpfr_clear(x);
}