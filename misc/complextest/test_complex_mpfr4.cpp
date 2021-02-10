#include <stdio.h>
#include <mpcomplex.h>
#include <complex>
#include <iostream>

void test_constructor()
{
   mpfr::mpcomplex a;
   mpfr::mpreal ra, rb, rc;
   ra = "1.0"; rb = "2.0"; rc = "3.0";
   mpfr::mpcomplex b(ra,rb);
   mpfr::mpcomplex c(rc);
   printf("constructor works\n");
}

void test_extraction()
{
   mpfr::mpcomplex a;
   mpfr::mpreal ra, rb;
   ra = "1.0"; rb = "2.0";
   mpfr::mpcomplex b(ra,rb);
   std::cout << b.real() << "\n";
   std::cout << b.imag() << "\n";
   b.real() = "3.14";
   b.imag() = "2.71";
   std::cout << b.real() << "\n";
   std::cout << b.imag() << "\n";
}

void test_comparison()
{
   mpfr::mpcomplex a (1.0, 2.0);
   mpfr::mpcomplex b (1.0, 2.0);
   mpfr::mpcomplex c (3.0, 4.0);
   mpfr::mpreal    d = 3.0;
   mpfr::mpcomplex e;
   e = 3.0;

  if (a==b) {printf("comparision1 ok\n");}
  if (b==a) {printf("comparision1 ok\n");}
  if (a!=b) {printf("comparision1 ng\n");}
  if (b!=a) {printf("comparision1 ng\n");}

  if (a==c) {printf("comparision2 ng\n");}
  if (c==a) {printf("comparision2 ng\n");}
  if (a!=c) {printf("comparision2 ok\n");}
  if (c!=a) {printf("comparision2 ok\n");}

  if (d==e) {printf("comparision3 ok\n");}
  if (e==d) {printf("comparision3 ok\n");}
  if (d!=e) {printf("comparision3 ng\n");}
  if (e!=d) {printf("comparision3 ng\n");}
}


int main(int argc, char *argv[])
{
  test_constructor();
  test_extraction();
  test_comparison();
  return 0;
}
