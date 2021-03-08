#include <mplapack.h>
#include <mplapack_debug.h>
#include <fstream>
#include <iostream>
#include <string>

extern "C" int dlasq1_ (int *n, double *d, double *e, double *work, int *info);

//Matlab/Octave format
void printvec (int N, mpreal * v) {
  mpreal mtmp;
  printf ("[ ");
  for (int i = 0; i < N; i++){
      mtmp = v[i];
      mpfr_printf ("%20.16Re", mpfr_ptr (mtmp));
      if (i < N - 1) printf (", ");
  }
  printf ("]; ");
}

void printvec (int N, double *v) {
  double mtmp;
  printf ("[ ");
  for (int i = 0; i < N; i++) {
      mtmp = v[i];
      printf ("%20.16e", mtmp);
      if (i < N - 1) printf (", ");
  }
  printf ("]; ");
}

int main ()
{
  std::ifstream ifs ("test_Rlasq1.txt");
  std::string str;
  while (getline (ifs, str)) {
      if (str.empty()) exit(1);
      int n = std::stoi(str);
      printf("#length of n is %d \n",n);
      mplapackint info;
      mpreal *d = new mpreal[n];
      mpreal *e = new mpreal[n - 1];
      mpreal *work = new mpreal[4 * n];
      int default_prec = 256;
      mpfr_set_default_prec (default_prec);

      int infod;
      double *dd = new double[n];
      double *ed = new double[n - 1];
      double *workd = new double[4 * n];
//setting d, e vector
      int i=0;
      while (i<n) {
          getline (ifs, str);
   	  d[i] = dd[i] = std::stod(str);
          i++;
      }
      i=0;
      while (i<n-1) {
          getline (ifs, str);
	  e[i] = ed[i] = std::stod(str);
          i++;
      }
      printf ("d ="); printvec (n, d);
      printf ("e ="); printvec (n - 1, e); printf ("\n");
//get singular values
      Rlasq1 (n, d, e, work, &info);
      dlasq1_ (&n, dd, ed, workd, &infod);
//print out some results.
      printf ("#singularvalues by mplapack\n");
      printvec (n, d); printf ("\n");
      printf ("#singularvalues by lapack \n");
      printvec (n, dd);  printf ("\n");

      double diff = infnorm (dd, d, n, 1);
      if (diff > EPSILON) { printf ("error: %lf", diff); printf ("\n"); exit(1); }
      else { printf ("ok: %lf", diff); printf ("\n"); }
//
      delete[]workd;
      delete[]dd;
      delete[]ed;
      delete[]work;
      delete[]d;
      delete[]e;
    }
}
