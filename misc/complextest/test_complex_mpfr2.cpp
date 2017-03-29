#include <stdio.h>
#include <mpcomplex.h>
#include <complex>
#include <iostream>

using namespace mpfr;
using namespace std;

void test_multiplication1()
{
  cout << "testing multiplication\n";
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  a *= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "multiplication works\n";
}

void test_multiplication2()
{
  cout << "testing multiplication\n";
  mpcomplex a(1.0, 2.0);
  mpc_t b;  mpc_init2 (b, 200);  mpc_set_d_d (b, 3.0, 4.0, GMP_RNDD);
  a *= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "multiplication works\n";
}

void test_multiplication3()
{
  cout << "testing multiplication\n";
  mpcomplex a(1.0, 2.0);
  complex<double> b (3.0, 4.0);
  a *= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "multiplication works\n";
}

void test_multiplication4()
{
  cout << "testing multiplication 4\n";
  mpcomplex a(-1.0, -2.0); 
  mpreal b = 1.1;
  a *= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "multiplication works\n";
}

void test_multiplication5()
{
  cout << "testing multiplication 5\n";
  mpcomplex a(-1.0, -2.0);
  mpfr_t s;  mpfr_init2 (s, 200); mpfr_set_d (s, 2.1, GMP_RNDD);
  a *= s;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "multiplication works\n";
}

void test_multiplication6()
{
  cout << "testing multiplication 6\n";
  mpcomplex a(-1.0, -2.0);
  double b = 1.111;
  a *= b;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "multiplication works\n";
}

void test_division1()
{
  cout << "testing division\n";
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  complex<double> ad(1.0, 2.0);
  complex<double> bd(3.0, 4.0);
  a /= b;
  ad /= bd;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";

  cout << ad.real() << "\n";
  cout << ad.imag() << "\n";
  cout << "division works\n";
}

void test_division2()
{
  cout << "testing division\n";
  mpcomplex a(1.0, 2.0);
  mpc_t b;  mpc_init2 (b, 200);  mpc_set_d_d (b, 3.0, 4.0, GMP_RNDD);
  a /= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "division works\n";
}

void test_division3()
{
  cout << "testing division\n";
  mpcomplex a(1.0, 2.0);
  complex<double> b (3.0, 4.0);
  a /= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "division works\n";
}

void test_division4()
{
  cout << "testing division 4\n";
  mpcomplex a(1.0, 2.0); 
  mpreal b = 0.5;
  a /= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "division works\n";
}

void test_division5()
{
  cout << "testing division 5\n";
  mpcomplex a(1.0, 2.0);
  mpfr_t s;  mpfr_init2 (s, 200); mpfr_set_d (s, 0.5, GMP_RNDD);
  a /= s;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "division works\n";
}

void test_equality()
{
  mpcomplex a(1.0, 2.0);
  mpcomplex b(1.0, 2.0);
  mpcomplex c(-1.0, 2.0);
  if (a==b) { cout << "a == b is ok\n"; }
  else { cout << "a == b is ng\n"; }
  if (a==c) { cout << "a == c is ng\n"; }
  else { cout << "a == c is ok\n"; }
}

void test_non_equality()
{
  mpcomplex a(1.0, 2.0);
  mpcomplex b(1.0, 2.0);
  mpcomplex c(-1.0, 2.0);
  if (a!=b) { cout << "a != b is ng\n"; }
  else { cout << "a != b is ok\n"; }
  if (a!=c) { cout << "a != c is ok\n"; }
  else { cout << "a != c is ng\n"; }
}

void test_division6()
{
  cout << "testing division 6\n";
  mpcomplex a(1.0, 2.0);
  double b = 0.5;
  a /= b;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "division works\n";
}

void test_cast2complexdouble()
{
  cout << "testing cast2complexdouble\n";
  mpcomplex a(1.0, 2.0);
  complex<double> b (0.0, 0.0);
  b = complex<double>(a);
  cout << b.real() << "\n";
  cout << b.imag() << "\n";
  cout << "cast works\n";
}

void test_pow()
{
  mpcomplex a (1.0, 3.0);
  mpcomplex b (3.0, 4.0);
  mpcomplex c (0.0, 0.0);

  complex<double> ad (1.0, 3.0);
  complex<double> bd (3.0, 4.0);
  complex<double> cd (0.0, 0.0);
  c = pow(a,b);
  cd = pow(a,b);

  cout << a.real() << "\n";
  cout << a.imag() << "\n";

  cout << c.real() << "\n";
  cout << c.imag() << "\n";

  cout << cd.real() << "\n";
  cout << cd.imag() << "\n";

  cout << "pow works\n";
}

void test_sqrt()
{
  mpcomplex a (1.0, 3.0);
  mpcomplex b;

  complex<double> ad (1.0, 3.0);
  complex<double> bd;

  b = sqrt(a);
  bd = sqrt(ad);

  cout << b.real() << "\n";
  cout << b.imag() << "\n";

  cout << bd.real() << "\n";
  cout << bd.imag() << "\n";

  cout << "sqrt works\n";
}

void test_abs()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << abs(a) << "\n";
  cout << abs(ad) << "\n";

  cout << "abs works\n";
}

void test_arg()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << arg(a) << "\n";
  cout << arg(ad) << "\n";

  cout << "arg works\n";
}

void test_proj()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << proj(a).real() << "\n";
  cout << proj(a).imag() << "\n";
  cout << proj(ad).real() << "\n";
  cout << proj(ad).imag() << "\n";

  cout << "proj works\n";
}

void test_conjg()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << conjg(a).real() << "\n";
  cout << conjg(a).imag() << "\n";
  cout << conjg(ad).real() << "\n";
  cout << conjg(ad).imag() << "\n";

  cout << "conjg works\n";
}

void test_norm()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << norm(a) << "\n";
  cout << norm(ad) << "\n";

  cout << "norm works\n";
}

void test_exp()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << exp(a).real() << "\n";
  cout << exp(a).imag() << "\n";

  cout << exp(ad).real() << "\n";
  cout << exp(ad).imag() << "\n";

  cout << "exp works\n";
}

void test_log()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << log(a).real() << "\n";
  cout << log(a).imag() << "\n";

  cout << log(ad).real() << "\n";
  cout << log(ad).imag() << "\n";

  cout << "log works\n";
}

void test_sin()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << sin(a).real() << "\n";
  cout << sin(a).imag() << "\n";

  cout << sin(ad).real() << "\n";
  cout << sin(ad).imag() << "\n";

  cout << "sin works\n";
}

void test_cos()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << cos(a).real() << "\n";
  cout << cos(a).imag() << "\n";

  cout << cos(ad).real() << "\n";
  cout << cos(ad).imag() << "\n";

  cout << "cos works\n";
}
void test_tan()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << tan(a).real() << "\n";
  cout << tan(a).imag() << "\n";

  cout << tan(ad).real() << "\n";
  cout << tan(ad).imag() << "\n";

  cout << "tan works\n";
}

void test_sinh()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << sinh(a).real() << "\n";
  cout << sinh(a).imag() << "\n";

  cout << sinh(ad).real() << "\n";
  cout << sinh(ad).imag() << "\n";

  cout << "sinh works\n";
}

void test_cosh()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << cosh(a).real() << "\n";
  cout << cosh(a).imag() << "\n";

  cout << cosh(ad).real() << "\n";
  cout << cosh(ad).imag() << "\n";

  cout << "cosh works\n";
}
void test_tanh()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << tanh(a).real() << "\n";
  cout << tanh(a).imag() << "\n";

  cout << tanh(ad).real() << "\n";
  cout << tanh(ad).imag() << "\n";

  cout << "tanh works\n";
}


void test_asin()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << asin(a).real() << "\n";
  cout << asin(a).imag() << "\n";

  cout << asin(ad).real() << "\n";
  cout << asin(ad).imag() << "\n";

  cout << "asin works\n";
}

void test_acos()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << acos(a).real() << "\n";
  cout << acos(a).imag() << "\n";

  cout << acos(ad).real() << "\n";
  cout << acos(ad).imag() << "\n";

  cout << "acos works\n";
}
void test_atan()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << atan(a).real() << "\n";
  cout << atan(a).imag() << "\n";

  cout << atan(ad).real() << "\n";
  cout << atan(ad).imag() << "\n";

  cout << "atan works\n";
}

void test_asinh()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << asinh(a).real() << "\n";
  cout << asinh(a).imag() << "\n";

  cout << asinh(ad).real() << "\n";
  cout << asinh(ad).imag() << "\n";

  cout << "asinh works\n";
}

void test_acosh()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << acosh(a).real() << "\n";
  cout << acosh(a).imag() << "\n";

  cout << acosh(ad).real() << "\n";
  cout << acosh(ad).imag() << "\n";

  cout << "acosh works\n";
}
void test_atanh()
{
  mpcomplex a (1.0, 3.0);
  complex<double> ad (1.0, 3.0);

  cout << atanh(a).real() << "\n";
  cout << atanh(a).imag() << "\n";

  cout << atanh(ad).real() << "\n";
  cout << atanh(ad).imag() << "\n";

  cout << "atanh works\n";
}


int main(int argc, char *argv[])
{
  test_pow();
  test_abs();
  test_sqrt();
  test_conjg();
  test_norm();
  test_arg();
  test_exp();
  test_log();

  test_sin();
  test_cos();
  test_tan();
  test_sinh();
  test_cosh();
  test_tanh();

  test_asin();
  test_acos();
  test_atan();
  test_asinh();
  test_acosh();
  test_atanh();
/*
  test_division1();
  test_division2();
  test_division3();
  test_division4();
  test_division5();
  test_division6();

  test_multiplication1();
  test_multiplication2();
  test_multiplication3();
  test_multiplication4();
  test_multiplication5();
  test_multiplication6();
  test_equality();
  test_non_equality();
  test_cast2complexdouble();
*/
  return 0;
}
