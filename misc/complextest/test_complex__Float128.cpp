#define ___MPLAPACK_BUILD_WITH_BINARY128___

#include <stdio.h>
#include <mpcomplex.h>
#include <complex>
#include <iostream>


using namespace mpfr;
using namespace std;

void test_assignment10()
{
  cout << "testing assignment\n";
  mpreal  b, c;
  _Float128 bb, cc;
  b = 5.0;
  bb = 2.0;
  cc = 3.0;
  c = b - bb ;
  cout << c << "\n";  
}

/*
void test_subtraction1()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  mpcomplex c(0.0, 0.0);
  c = a - b;
  cout << c.real() <<  " + " << c.imag() << "i\n";

  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 90;
  mpcomplex d(1.0, 2.0);

  mpcomplex::default_real_prec = 20;
  mpcomplex::default_imag_prec = 10;
  mpcomplex e(3.0, 4.0);

  mpcomplex f = d - e;
  cout << f.real() <<  " + " << f.imag() << "i\n";
   
  cout << f.get_prec_re() << "\n";
  cout << f.get_prec_im() << "\n";

  cout << "subtraction1 works\n";
}

void test_subtraction2()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  std::complex<double> b(3.0, 4.0);
  mpcomplex c = a - b;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction2 works\n";
}

void test_subtraction3()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  const char *s = "(3.1415 1.1123)";
  mpcomplex c = a - s;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction3 works\n";
}

void test_subtraction4()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 100;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  std::complex<double> b(3.9, 4.1);
  mpcomplex c = b - a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction4 works\n";
}

void test_subtraction5()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 100;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  const char *s = "(3.1415926 2.71727)";
  mpcomplex c =  s - a; 
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction5 works\n";
}

void test_subtraction6()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 100;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  std::complex<double> b(4.0, 4.0);
  mpcomplex c = b - a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction6 works\n";
}

void test_subtraction7()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  mpreal::default_prec = 30;
  mpreal b (100);
  mpcomplex c = a - b;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  mpcomplex::default_real_prec = 50;
  mpcomplex::default_imag_prec = 50;
  mpreal::default_prec = 40;
  mpcomplex u(1.0, 2.0);

  mpreal v (100);
  mpcomplex w = u - v;

  cout << w.real() <<  " + " << w.imag() << "i\n";
  cout << w.get_prec_re() << "\n";
  cout << w.get_prec_im() << "\n";

  mpcomplex::default_real_prec = 40;
  mpcomplex::default_imag_prec = 40;
  mpreal::default_prec = 40;
  mpcomplex p(1.0, 2.0);

  mpreal q (100);
  mpcomplex r = p - q;

  cout << r.real() <<  " + " << r.imag() << "i\n";
  cout << r.get_prec_re() << "\n";
  cout << r.get_prec_im() << "\n";

  cout << "subtraction7 works\n";
}

void test_subtraction8()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  double b=3;
  mpcomplex c = a - b;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction8 works\n";
}
void test_subtraction9()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  int b=3;
  mpcomplex c = a - b;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction9 works\n";
}

void test_subtraction10()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  double b=3;
  mpcomplex c = b - a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction10 works\n";
}

void test_subtraction11()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  int b=3;
  mpcomplex c = b - a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "subtraction10 works\n";
}
*/

int main(int argc, char *argv[])
{
  test_assignment10();
/*
  test_subtraction1();
  test_subtraction2();
  test_subtraction3();
  test_subtraction4();
  test_subtraction5();
  test_subtraction6();
  test_subtraction7();
  test_subtraction8();
  test_subtraction9();
  test_subtraction10();
  test_subtraction11();
*/
  return 0;
}
