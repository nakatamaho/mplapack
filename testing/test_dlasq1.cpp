#include <iostream>
#include <stdio.h>

extern "C" int dlasq1_(int *n, double *d, double *e, double *work, int *info);

//Matlab/Octave format
void printvec(int N, double *v) {
  double mtmp;
  printf("[ ");
  for (int i = 0; i < N; i++) {
    mtmp = v[i];
    printf("%5.2e", mtmp);
    if (i < N-1 ) printf (", ");
  } 
  printf("]; ");
}
int main()
{
  int n = 3;
  int info;
  double *d = new double[n];
  double *e = new double[n-1];
  double *work = new double[4*n];
//setting d, e vector
  d[0]=1;d[1]=2;d[2]=3;
  e[0]=1;e[1]=2;
  printf("d ="); printvec(n, d);
  printf("e ="); printvec(n-1, e);
  printf("\n");
//get singularvalue
  dlasq1_(&n, d, e, work, &info);
//print out some results.
  printf("#singularvalues \n"); printf("w =");
  printvec(n, d); printf("\n");
  delete[]work;
  delete[]d;
  delete[]e;
}
