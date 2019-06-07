#include <mblas_mpfr.h>
#include <mlapack_mpfr.h>

//Matlab/Octave format
void printvec(int N, mpreal *v) {
  mpreal mtmp;
  printf("[ ");
  for (int i = 0; i < N; i++) {
    mtmp = v[i];
    mpfr_printf("%8.6Re", mpfr_ptr(mtmp));
    if (i < N-1 ) printf (", ");
  } 
  printf("]; ");
}
int main()
{
  int n = 3;
  mpackint info;
  mpreal *d = new mpreal[n];
  mpreal *e = new mpreal[n-1];
  mpreal *work = new mpreal[4*n];
  int default_prec = 256;
  mpfr_set_default_prec(default_prec);

//setting d, e vector
  d[0]=1;d[1]=2;d[2]=3;
  e[0]=1;e[1]=2;
  printf("d ="); printvec(n, d);
  printf("e ="); printvec(n-1, e);
  printf("\n");
//get singularvalue
  Rlasq1(n, d, e, work, &info);
//print out some results.
  printf("#singularvalues \n"); printf("w =");
  printvec(n, d); printf("\n");
  delete[]work;
  delete[]d;
  delete[]e;
}
