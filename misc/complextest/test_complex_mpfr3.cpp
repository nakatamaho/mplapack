#include <stdio.h>
#include <mpcomplex.h>
#include <complex>
#include <iostream>

using namespace mpfr;
using namespace std;

void test_addition1()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  mpcomplex c(0.0, 0.0);
  c = a + b;
  cout << c.real() <<  " + " << c.imag() << "i\n";

  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 90;
  mpcomplex d(1.0, 2.0);

  mpcomplex::default_real_prec = 20;
  mpcomplex::default_imag_prec = 10;
  mpcomplex e(3.0, 4.0);

  mpcomplex f = d + e;
  cout << f.real() <<  " + " << f.imag() << "i\n";
   
  cout << f.get_prec_re() << "\n";
  cout << f.get_prec_im() << "\n";

  cout << "addition1 works\n";
}

void test_addition2()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  std::complex<double> b(3.0, 4.0);
  mpcomplex c = a + b;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition2 works\n";
}

void test_addition3()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  const char *s = "(3.1415 1.1123)";
  mpcomplex c = a + s;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition3 works\n";
}

void test_addition4()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 100;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  std::complex<double> b(3.9, 4.1);
  mpcomplex c = b + a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition4 works\n";
}

void test_addition5()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  const char *s = "(3.1415 1.1123)";
  mpcomplex c =  s + a; 
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition5 works\n";
}

void test_addition6()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 100;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  std::complex<double> b(4.0, 4.0);
  mpcomplex c = b + a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition6 works\n";
}

void test_addition7()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  mpreal::default_prec = 30;
  mpreal b (100);
  mpcomplex c = b + a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition6 works\n";
}

void test_addition8()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  double b=3;
  mpcomplex c = b + a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition8 works\n";
}
void test_addition9()
{
  cout << "testing addition\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 100;
  mpcomplex a(1.0, 2.0);
  int b=3;
  mpcomplex c = b + a;
  
  cout << c.real() <<  " + " << c.imag() << "i\n";
  cout << c.get_prec_re() << "\n";
  cout << c.get_prec_im() << "\n";

  cout << "addition9 works\n";
}


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

  cout << "subtraction works\n";
}

void test_subtraction2()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  mpreal b(3.0, 4.0);
  mpcomplex c(0.0, -2.0);
  c = a - b;
  cout << "subtraction works\n";
}


void test_multiplication1()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  mpcomplex c(0.0, 0.0);
  c = a * b;
  cout << c.real() <<  " + " << c.imag() << "i\n";

  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 90;
  mpcomplex d(1.0, 2.0);

  mpcomplex::default_real_prec = 20;
  mpcomplex::default_imag_prec = 10;
  mpcomplex e(3.0, 4.0);

  mpcomplex f = d * e;
  cout << f.real() <<  " + " << f.imag() << "i\n";
   
  cout << f.get_prec_re() << "\n";
  cout << f.get_prec_im() << "\n";

  cout << "multiplication works\n";
}


void test_division1()
{
  cout << "testing subtraction\n";
  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 10;
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  mpcomplex c(0.0, 0.0);
  c = a / b;
  cout.precision(30);
  cout << c.real() <<  " + " << c.imag() << "i\n";

  mpcomplex::default_real_prec = 10;
  mpcomplex::default_imag_prec = 90;
  mpcomplex d(1.0, 2.0);

  mpcomplex::default_real_prec = 20;
  mpcomplex::default_imag_prec = 10;
  mpcomplex e(3.0, 4.0);

  mpcomplex f = d / e; cout.precision(90);
  cout << f.real() <<  " + " << f.imag() << "i\n";
   
  cout << f.get_prec_re() << "\n";
  cout << f.get_prec_im() << "\n";

  cout << "division works\n";
}



void test_get_prec()
{
  cout << "testing get_prec\n";
  mpcomplex a(1.0, 2.0);
  cout << a.get_prec() << "\n";
  cout << "get_prec works\n";
}

void test_get_prec2()
{
  mpcomplex::default_real_prec = 102;
  mpcomplex::default_imag_prec = 204;
  cout << "testing get_prec\n";
  mpcomplex a(1.0, 2.0);
  cout << a.get_prec_re() << "\n";
  cout << a.get_prec_im() << "\n";
  cout << "get_prec_re and im works\n";
}

void test_set_prec()
{
  mpcomplex::default_real_prec = 102;
  mpcomplex::default_imag_prec = 204;
  cout << "testing set_prec\n";
  mpcomplex a(1.0, 2.0);
  cout << a.get_prec_re() << "\n";
  cout << a.get_prec_im() << "\n";
  a.set_prec(309);
  cout << a.get_prec_re() << "\n";
  cout << a.get_prec_im() << "\n";
  cout << a.real() << " + i" << a.imag() << "\n";

  a.set_prec2(10,30);
  cout << a.get_prec_re() << "\n";
  cout << a.get_prec_im() << "\n";
  cout << a.real() << " + i" << a.imag() << "\n";

  cout << "set_prec works\n";
}

int main(int argc, char *argv[])
{
  test_get_prec();
  test_get_prec2();
  test_set_prec();

  test_addition1();
  test_subtraction1();
  test_subtraction2();
  test_multiplication1();
  test_division1();

  test_addition2();
  test_addition3();
  test_addition4();
  test_addition5();
  test_addition6();
  test_addition7();
  test_addition8();
  test_addition9();
  return 0;
}
