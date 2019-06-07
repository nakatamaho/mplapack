#include <mlapack.h>
#include <mpack_debug.h>


extern "C" int dlasq1_(int *n, double *d, double *e, double *work, int *info);

//Matlab/Octave format
void printvec(int N, mpreal *v) {
  mpreal mtmp;
  printf("[ ");
  for (int i = 0; i < N; i++) {
    mtmp = v[i];
    mpfr_printf("%20.16Re", mpfr_ptr(mtmp));
    if (i < N-1 ) printf (", ");
  } 
  printf("]; ");
}

void printvec(int N, double *v) {
  double mtmp;
  printf("[ ");
  for (int i = 0; i < N; i++) {
    mtmp = v[i];
    printf("%20.16e", mtmp);
    if (i < N-1 ) printf (", ");
  } 
  printf("]; ");
}

int main()
{
  int n = 2;

  mpackint info;
  mpreal *d = new mpreal[n];
  mpreal *e = new mpreal[n-1];
  mpreal *work = new mpreal[4*n];
  int default_prec = 256;
  mpfr_set_default_prec(default_prec);

  int infod;
  double *dd = new double[n];
  double *ed = new double[n-1];
  double *workd = new double[4*n];

//setting d, e vector
  d[0]=dd[0]=0; d[1]=dd[1]=2; //d[2]=dd[2]=3;
  e[0]=ed[0]=0; // e[1]=ed[1]=2;

  printf("d ="); printvec(n, d);
  printf("e ="); printvec(n-1, e);
  printf("\n");
//get singularvalue
  Rlasq1(n, d, e, work, &info);
  dlasq1_(&n, dd, ed, workd, &infod);

//print out some results.
  printf("#singularvalues by mplapack\n");
  printvec(n, d); printf("\n");
  printf("#singularvalues by lapack \n");
  printvec(n, dd); printf("\n");

  double diff = infnorm(dd, d, n, 1);
  if (diff > EPSILON) {
      printf("error: %lf", diff); printf("\n");
  } else {
      printf("ok: %lf", diff); printf("\n");
  }

  delete[]workd;
  delete[]dd;
  delete[]ed;
  delete[]work;
  delete[]d;
  delete[]e;
}
