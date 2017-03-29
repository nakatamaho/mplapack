#include <stdio.h>
#include <mpcomplex.h>
#include <complex>
#include <iostream>

using namespace mpfr;
using namespace std;

void test_assignment10()
{
  cout << "testing assignment\n";
  mpcomplex a;
  mpreal  b = 1.0, c = -1.0 ;
  a.real() = b;
  a.real() = c;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "assignment works\n";
}

void test_constructor1()
{
  cout << "testing constructor mpcomplex(mpc_t)\n";
  mpc_t s;
  mpc_init2 (s, 200);
  mpc_set_d_d (s, 2.1, 1.3, GMP_RNDD);
  mpcomplex a(s);

  cout.precision(128);
  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(mpc_t) works\n";
} 

void test_constructor2()
{
  cout << "testing constructor mpcomplex(mpfr_t,  mpfr_t)\n";
  mpfr_t s, t;
  mpfr_init2 (s, 200);
  mpfr_set_d (s, 2.1, GMP_RNDD);

  mpfr_init2 (t, 200);
  mpfr_set_d (t, 1.11, GMP_RNDD);
  mpcomplex a(s,t);
  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(mpfr_t,  mpfr_t) works\n";
} 

void test_constructor3()
{
  cout << "testing constructor mpcomplex(mpf_t,  mpf_t)\n";
  mpf_t s, t, u;
  mpf_init2 (s, 200);
  mpf_set_d (s, 2.1);

  mpf_init2 (t, 200);
  mpf_set_d (t, 1.11);
  mpcomplex a(s,t);
  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(mpfr_t,  mpfr_t) works\n";
} 

void test_constructor4()
{
  cout << "testing constructor mpcomplex(char *)\n";
  mpcomplex a ("(3.14159 -7.01)");
  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(char*) works\n";
}

void test_constructor5()
{
  cout << "testing constructor mpcomplex(&mpcomplex)\n";
  mpcomplex a ("(3.14159 -7.01)");
  mpcomplex b (a);
  cout << b.real().to_string() << "\n";
  cout << b.imag().to_string() << "\n";
  cout << "constructor mpcomplex(&mpcomplex) is ok \n";
}

void test_constructor6()
{
  cout << "testing constructor mpcomplex(complex<double>)\n";
  complex <double> p (100.0,200.0);
  mpcomplex a(p);

  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(complex<double>) works \n";
}

void test_constructor7()
{
  cout << "testing constructor mpcomplex(mpreal, mpreal)\n";
   mpreal ra, rb;
   ra = "3.1415926"; 
   rb = "1.7320508";
   mpcomplex a(ra,rb);

  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(mpreal, mpreal) works \n";
}

void test_constructor8()
{
  cout << "testing constructor mpcomplex(double, double)\n";
   double ra = 1.100;
   double rb = 2.111;
   mpcomplex a(ra, rb);

  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(double, double) works \n";
}

void test_constructor9()
{
  cout << "testing constructor mpcomplex(char *, char *)\n";
  mpcomplex a ("3.14159", "-7.01");
  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(char*, char *) works\n";
}

void test_constructor10()
{
  cout << "testing constructor mpcomplex(mpreal)\n";
   mpreal ra, rb;
   ra = "3.1415926"; 
   mpcomplex a(ra);

  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor10 mpcomplex(mpreal) works \n";
}

void test_constructor11()
{
  cout << "testing constructor mpcomplex(mpfr_t)\n";
  mpfr_t s;
  mpfr_init2 (s, 200);
  mpfr_set_d (s, 2.1, GMP_RNDD);

  mpcomplex a(s);

  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(mpfr_t) works \n";
}

void test_constructor12()
{
  cout << "testing constructor mpcomplex(mpf_t)\n";
  mpf_t s;
  mpf_init2 (s, 200);
  mpf_set_d (s, 2.10000);

  mpcomplex a(s);

  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(mpf_t) works \n";
}

void test_constructor13()
{
  cout << "testing constructor mpcomplex(double)\n";
  double p = 1.001;
  mpcomplex a(p);

  cout.precision(128);
  cout << a.real().to_string() << "\n";
  cout << a.imag().to_string() << "\n";
  cout << "constructor mpcomplex(double) works \n";
}

void test_extraction()
{
  cout << "testing extraction\n";
  mpcomplex a ("(3.14159 -7.0199998)");
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "extraction works\n";
}

void test_assignment1()
{
  cout << "testing assignment\n";
  mpcomplex a ("(3.14159 -7.0199998)");
  mpcomplex b;
  b = a;
  cout.precision(100);
  cout << b.real() << "\n";
  cout << b.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment2()
{
  cout << "testing assignment\n";
  mpc_t s;  mpc_init2 (s, 200);  mpc_set_d_d (s, 2.1, 1.3, GMP_RNDD);
  mpcomplex b;
  b = s;
  cout.precision(100);
  cout << b.real() << "\n";
  cout << b.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment3()
{
  cout << "testing assignment\n";
  complex<double> s(1.0,-3.11);
  mpcomplex b;
  b = s;
  cout.precision(100);
  cout << b.real() << "\n";
  cout << b.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment4()
{
  cout << "testing assignment\n";
  const char *s = "(3.1415 1.1123)";
  mpcomplex b;
  b = s;
  cout.precision(100);
  cout << b.real() << "\n";
  cout << b.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment5()
{
  cout << "testing assignment\n";
  mpcomplex a ("(3.14159 -7.0199998)");
  mpreal b;
  b = 1.0;
  a = b;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment6()
{
  cout << "testing assignment\n";
  mpcomplex a ("(3.14159 -7.0199998)");
  mpfr_t s;  mpfr_init2 (s, 200);  mpfr_set_d (s, 2.1, GMP_RNDD);
  a = s;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment7()
{
  cout << "testing assignment\n";
  mpcomplex a ("(3.14159 -7.0199998)");
  double b = 1.0;
  a = b;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "assignment works\n";
}

void test_assignment8()
{
  cout << "testing assignment\n";
  mpcomplex a (3.14159, -7.0199998);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "assignment8 works\n";
}


void test_assignment9()
{
  cout << "testing assignment\n";
  mpcomplex a;
  a = mpcomplex(3.14159, -7.0199998);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "assignment9 works\n";
}


void test_addition1()
{
  cout << "testing addition\n";
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  a += b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition2()
{
  cout << "testing addition\n";
  mpcomplex a(1.0, 2.0);
  mpc_t b;  mpc_init2 (b, 200);  mpc_set_d_d (b, 3.0, 4.0, GMP_RNDD);
  a += b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition3()
{
  cout << "testing addition\n";
  mpcomplex a(1.0, 2.0);
  complex<double> b (3.0, 4.0);
  a += b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition4()
{
  cout << "testing addition 4\n";
  mpcomplex a(-1.0, -2.0); 
  mpreal b = 1.1;
  a += b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition5()
{
  cout << "testing addition 5\n";
  mpcomplex a(-1.0, -2.0);
  mpfr_t s;  mpfr_init2 (s, 200); mpfr_set_d (s, 2.1, GMP_RNDD);
  a += s;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition6()
{
  cout << "testing addition 6\n";
  mpcomplex a(-1.0, -2.0);
  a = +a;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition7()
{
  cout << "testing addition 7\n";
  mpcomplex a(-1.0, -2.0);
  double b = 1.111;
  a += b;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition8()
{
  cout << "testing addition 8\n";
  mpcomplex a(-1.0, -2.0);
  a++;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_addition9()
{
  cout << "testing addition 9\n";
  mpcomplex a(-1.0, -2.0);
  ++a;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "addition works\n";
}

void test_subtraction1()
{
  cout << "testing subtraction\n";
  mpcomplex a(1.0, 2.0);
  mpcomplex b(3.0, 4.0);
  a -= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction2()
{
  cout << "testing subtraction\n";
  mpcomplex a(1.0, 2.0);
  mpc_t b;  mpc_init2 (b, 200);  mpc_set_d_d (b, 3.0, 4.0, GMP_RNDD);
  a -= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction3()
{
  cout << "testing subtraction\n";
  mpcomplex a(1.0, 2.0);
  complex<double> b (3.0, 4.0);
  a -= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction4()
{
  cout << "testing subtraction 4\n";
  mpcomplex a(-1.0, -2.0); 
  mpreal b = 1.1;
  a -= b;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction5()
{
  cout << "testing subtraction 5\n";
  mpcomplex a(-1.0, -2.0);
  mpfr_t s;  mpfr_init2 (s, 200); mpfr_set_d (s, 2.1, GMP_RNDD);
  a -= s;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction6()
{
  cout << "testing subtraction 6\n";
  mpcomplex a(-1.0, -2.0);
  a = -a;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction7()
{
  cout << "testing subtraction 7\n";
  mpcomplex a(-1.0, -2.0);
  double b = 1.111;
  a -= b;
  cout.precision(100);
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction8()
{
  cout << "testing subtraction 8\n";
  mpcomplex a(-1.0, -2.0);
  a--;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}

void test_subtraction9()
{
  cout << "testing subtraction 9\n";
  mpcomplex a(-1.0, -2.0);
  --a;
  cout << a.real() << "\n";
  cout << a.imag() << "\n";
  cout << "subtraction works\n";
}


int main(int argc, char *argv[])
{

  test_constructor1();
  test_constructor2();
  test_constructor3();
  test_constructor4();
  test_constructor5();
  test_constructor6();
  test_constructor7();
  test_constructor8();
  test_constructor9();
  test_constructor10();
  test_constructor11();
  test_constructor12();
  test_constructor13();
  test_extraction();
  test_assignment1();
  test_assignment2();
  test_assignment3();
  test_assignment4();
  test_assignment5();
  test_assignment6();
  test_assignment7();
  test_assignment8();
  test_assignment9();

  test_addition1();
  test_addition2();
  test_addition3();
  test_addition4();
  test_addition5();
  test_addition6();
  test_addition7();
  test_addition8();
  test_addition9();

  test_subtraction1();
  test_subtraction2();
  test_subtraction3();
  test_subtraction4();
  test_subtraction5();
  test_subtraction6();
  test_subtraction7();
  test_subtraction8();
  test_subtraction9();


  return 0;
}
