/*
 * Copyright (c) 2008-2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
#include <complex>
#include <mpblas.h>
#include <mplapack_debug.h>
#include <string.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MAX_ITER 10

void subst_test1()
{
  char buf1[BUFLEN], buf2[BUFLEN], buf3[BUFLEN];
  REAL tmp1;
  REAL_REF tmp2;
  printf("*** Substitution test 1 ***\n");
  strcpy (buf1, "-1.234567890123456789012345678901234567890123456789012345678901234567890e+01");

//tmp1 = buf1;
#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
  tmp1 = strtoflt128(buf1, NULL);
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
  sscanf(buf1, "%lf", &tmp1);
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
  sscanf(buf1, "%Lf", &tmp1);
#else
  tmp1 = buf1;
#endif

  tmp2 = tmp1;

  sprintnum(buf2, tmp1);
  sprintnum(buf3, tmp2);

#if defined VERBOSE_TEST
  printf("original  :%s\n",buf1);
  printf("mplib     :%s\n",buf2);
  printf("subst2refm:%s\n",buf3);
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
  if (strncmp(buf1, buf2, 19)==0 && strncmp(buf2, buf3, 19)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
  if (strncmp(buf1, buf2, 65)==0 && strncmp(buf2, buf3, BUFLEN)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_DD___
  if (strncmp(buf1, buf2, 34)==0 && strncmp(buf2, buf3, 34)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_QD___
  if (strncmp(buf1, buf2, 66)==0 && strncmp(buf2, buf3, 66)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
  if (strncmp(buf1, buf2, 19)==0 && strncmp(buf2, buf3, 19)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT128___
  if (strncmp(buf1, buf2, 37)==0 && strncmp(buf2, buf3, 37)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#endif
  printf("*** Substitution test 1 successful ***\n");
}

void subst_test2()
{
  char buf1[BUFLEN], buf2[BUFLEN], buf3[BUFLEN];
  REAL tmp1;
  REAL_REF tmp2;
  printf("*** Substitution test 2 ***\n");
  strcpy (buf1, "-1.234567890123456789012345678901234567890123456789012345678901234567890e+01");

//tmp2 = buf1;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
  sscanf(buf1, "%lf", &tmp2);
#else
  tmp2 = buf1;
#endif

//tmp1 = tmp2;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
  tmp1 = tmp2;
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
  tmp1 = cast2mpf_class(tmp2);
#elif defined ___MPLAPACK_BUILD_WITH_DD___
  tmp1 = cast2dd_real(tmp2);
#elif defined ___MPLAPACK_BUILD_WITH_QD___
  tmp1 = cast2qd_real(tmp2);
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
  tmp1 = cast2double(tmp2);
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
  tmp1 = cast2longdouble(tmp2);  
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT128___
  tmp1 = cast2_Float128(tmp2);
#endif

  sprintnum(buf2, tmp2);
  sprintnum(buf3, tmp1);
#if defined VERBOSE_TEST
  printf("original  :%s\n",buf1);
  printf("referencem:%s\n",buf2);
  printf("castedtomp:%s\n",buf3);
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
  if (strncmp(buf1, buf2, 19)==0 && strncmp(buf2, buf3, 19)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
  if (strncmp(buf1, buf2, 65)==0 && strncmp(buf2, buf3, BUFLEN)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_DD___
  if (strncmp(buf1, buf2, 34)==0 && strncmp(buf2, buf3, 34)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_QD___
  if (strncmp(buf1, buf2, 66)==0 && strncmp(buf2, buf3, 66)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
  if (strncmp(buf1, buf2, 19)==0 && strncmp(buf2, buf3, 19)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);}
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
  if (strncmp(buf1, buf2, 19)==0 && strncmp(buf2, buf3, 19)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);}
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT128___
  if (strncmp(buf1, buf2, 36)==0 && strncmp(buf2, buf3, 36)==0)  printf("ok!\n"); else {printf("failed!\n"); exit(1);} 
#endif
  printf("*** Substitution test 2 successful ***\n");
}

void mp_sub_test_complex()
{
    COMPLEX_REF Ctemp1r, Ctemp2r, dummy;
    COMPLEX Ctemp1;
    REAL_REF diff;

    set_random_number(Ctemp1r, Ctemp1);

    cout << "C1R  = "; printnum(Ctemp1r); cout << endl;
    cout << "C1   = "; printnum(Ctemp1);  cout << endl;

    Ctemp2r = Ctemp1 - Ctemp1r;
    diff = abs(Ctemp2r);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("error\n"); exit(1); }
}

void mp_sub_test_real()
{
    REAL_REF Rtemp1r, Rtemp2r, Rtemp3r;
    REAL Rtemp1;
    REAL_REF diff;

    printf("Substitution test\n");
    set_random_number(Rtemp1r, Rtemp1);
    printnum(Rtemp1r); printf("\n");
    printnum(Rtemp1);  printf("\n");
    Rtemp2r = Rtemp1 - Rtemp1r;
    printnum(Rtemp2r); printf("\n");
    Rtemp3r = Rtemp1r - Rtemp1;
    printnum(Rtemp3r); printf("\n");

    diff = abs(Rtemp2r);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("error\n"); exit(1); }

    diff = abs(Rtemp3r);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("error\n"); exit(1); }
    printf("Substitution test successful\n");

    printf("Subtraction test\n");
    Rtemp1r = 1.0;
    Rtemp1 = 2.0;
    Rtemp2r = Rtemp1r - Rtemp1;    
    printnum(Rtemp2r);  printf("\n");

    Rtemp2r = Rtemp1 - Rtemp1r;
    printnum(Rtemp2r);  printf("\n");
    printf("Subtraction test successful\n");
}

void addition_real_test()
{
    REAL_REF Rtemp1r, Rtemp2r, Rtemp3r;
    REAL Rtemp1, Rtemp2, Rtemp3;
    REAL_REF diff;

    printf("*** REAL Addition test ***\n");
    set_random_number(Rtemp1r, Rtemp1);
    set_random_number(Rtemp2r, Rtemp2);

    Rtemp3r = Rtemp1r + Rtemp2r;
    Rtemp3  = Rtemp1  + Rtemp2;

    printnum(Rtemp1r); printf("\n");
    printnum(Rtemp1); printf("\n");
    printf("+"); printf("\n");
    printnum(Rtemp2r); printf("\n");
    printnum(Rtemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Rtemp3r); printf("\n");
    printnum(Rtemp3); printf("\n");

    diff = abs(Rtemp3r - Rtemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** REAL Addition test failed ***\n"); exit(1); }
    printf("*** REAL Addition test successful ***\n");
}

void subtraction_real_test()
{
    REAL_REF Rtemp1r, Rtemp2r, Rtemp3r;
    REAL Rtemp1, Rtemp2, Rtemp3;
    REAL_REF diff;

    printf("*** REAL Subtraction test ***\n");
    set_random_number(Rtemp1r, Rtemp1);
    set_random_number(Rtemp2r, Rtemp2);

    Rtemp3r = Rtemp1r - Rtemp2r;
    Rtemp3  = Rtemp1  - Rtemp2;

    printnum(Rtemp1r); printf("\n");
    printnum(Rtemp1); printf("\n");
    printf("-"); printf("\n");
    printnum(Rtemp2r); printf("\n");
    printnum(Rtemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Rtemp3r); printf("\n");
    printnum(Rtemp3); printf("\n");

    diff = abs(Rtemp3r - Rtemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** REAL Subtraction test failed ***\n"); exit(1); }
    printf("*** REAL Subtraction test successful ***\n");
}

void multiplication_real_test()
{
    REAL_REF Rtemp1r, Rtemp2r, Rtemp3r;
    REAL Rtemp1, Rtemp2, Rtemp3;
    REAL_REF diff;

    printf("*** REAL Multiplication test ***\n");
    set_random_number(Rtemp1r, Rtemp1);
    set_random_number(Rtemp2r, Rtemp2);

    Rtemp3r = Rtemp1r * Rtemp2r;
    Rtemp3  = Rtemp1  * Rtemp2;

    printnum(Rtemp1r); printf("\n");
    printnum(Rtemp1); printf("\n");
    printf("*"); printf("\n");
    printnum(Rtemp2r); printf("\n");
    printnum(Rtemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Rtemp3r); printf("\n");
    printnum(Rtemp3); printf("\n");

    diff = abs(Rtemp3r - Rtemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** REAL Multiplication test failed ***\n"); exit(1); }
    printf("*** REAL Multiplication test successful ***\n");
}

void division_real_test()
{
    REAL_REF Rtemp1r, Rtemp2r, Rtemp3r;
    REAL Rtemp1, Rtemp2, Rtemp3;
    REAL_REF diff;

    printf("*** REAL Division test ***\n");
    set_random_number(Rtemp1r, Rtemp1);
    set_random_number(Rtemp2r, Rtemp2);

    Rtemp3r = Rtemp1r / Rtemp2r;
    Rtemp3  = Rtemp1  / Rtemp2;

    printnum(Rtemp1r); printf("\n");
    printnum(Rtemp1); printf("\n");
    printf("/"); printf("\n");
    printnum(Rtemp2r); printf("\n");
    printnum(Rtemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Rtemp3r); printf("\n");
    printnum(Rtemp3); printf("\n");

    diff = abs(Rtemp3r - Rtemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** REAL Division test failed ***\n"); exit(1); }
    printf("*** REAL Division test successful ***\n");
}

void addition_complex_test()
{
    COMPLEX_REF Ctemp1r, Ctemp2r, Ctemp3r;
    COMPLEX Ctemp1, Ctemp2, Ctemp3;
    COMPLEX_REF diff;

    printf("*** COMPLEX Addition test ***\n");
    set_random_number(Ctemp1r, Ctemp1);
    set_random_number(Ctemp2r, Ctemp2);

    Ctemp3r = Ctemp1r + Ctemp2r;
    Ctemp3  = Ctemp1  + Ctemp2;

    printnum(Ctemp1r); printf("\n");
    printnum(Ctemp1); printf("\n");
    printf("+"); printf("\n");
    printnum(Ctemp2r); printf("\n");
    printnum(Ctemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Ctemp3r); printf("\n");
    printnum(Ctemp3); printf("\n");

    diff = abs(Ctemp3r - Ctemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** COMPLEX Addition test failed ***\n"); exit(1); }
    printf("*** COMPLEX Addition test successful ***\n");
}

void subtraction_complex_test()
{
    COMPLEX_REF Ctemp1r, Ctemp2r, Ctemp3r;
    COMPLEX Ctemp1, Ctemp2, Ctemp3;
    COMPLEX_REF diff;

    printf("*** COMPLEX Subtraction test ***\n");
    set_random_number(Ctemp1r, Ctemp1);
    set_random_number(Ctemp2r, Ctemp2);

    Ctemp3r = Ctemp1r - Ctemp2r;
    Ctemp3  = Ctemp1  - Ctemp2;

    printnum(Ctemp1r); printf("\n");
    printnum(Ctemp1); printf("\n");
    printf("-"); printf("\n");
    printnum(Ctemp2r); printf("\n");
    printnum(Ctemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Ctemp3r); printf("\n");
    printnum(Ctemp3); printf("\n");

    diff = abs(Ctemp3r - Ctemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** COMPLEX Subtraction test failed ***\n"); exit(1); }
    printf("*** COMPLEX Subtraction test successful ***\n");
}

void multiplication_complex_test()
{
    COMPLEX_REF Ctemp1r, Ctemp2r, Ctemp3r;
    COMPLEX Ctemp1, Ctemp2, Ctemp3;
    COMPLEX_REF diff;

    printf("*** COMPLEX Multiplication test ***\n");
    set_random_number(Ctemp1r, Ctemp1);
    set_random_number(Ctemp2r, Ctemp2);

    Ctemp3r = Ctemp1r * Ctemp2r;
    Ctemp3  = Ctemp1  * Ctemp2;

    printnum(Ctemp1r); printf("\n");
    printnum(Ctemp1); printf("\n");
    printf("*"); printf("\n");
    printnum(Ctemp2r); printf("\n");
    printnum(Ctemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Ctemp3r); printf("\n");
    printnum(Ctemp3); printf("\n");

    diff = abs(Ctemp3r - Ctemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** COMPLEX Multiplication test failed ***\n"); exit(1); }
    printf("*** COMPLEX Multiplication test successful ***\n");
}

void division_complex_test()
{
    COMPLEX_REF Ctemp1r, Ctemp2r, Ctemp3r;
    COMPLEX Ctemp1, Ctemp2, Ctemp3;
    COMPLEX_REF diff;

    printf("*** COMPLEX Division test ***\n");
    set_random_number(Ctemp1r, Ctemp1);
    set_random_number(Ctemp2r, Ctemp2);

    Ctemp3r = Ctemp1r / Ctemp2r;
    Ctemp3  = Ctemp1  / Ctemp2;

    printnum(Ctemp1r); printf("\n");
    printnum(Ctemp1); printf("\n");
    printf("/"); printf("\n");
    printnum(Ctemp2r); printf("\n");
    printnum(Ctemp2); printf("\n");
    printf("="); printf("\n");
    printnum(Ctemp3r); printf("\n");
    printnum(Ctemp3); printf("\n");

    diff = abs(Ctemp3r - Ctemp3);
    cout << "diff = "; printnum(diff); cout << endl;
    if (abs(diff) > EPSILON) { printf("*** COMPLEX Division test failed ***\n"); exit(1); }
    printf("*** COMPLEX Division test successful ***\n");
}

int main(int argc, char *argv[])
{
    printf("*** Testing arithmetic start ***\n");

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    mpf_set_default_prec(___MPLAPACK_DEFAULT_PRECISION___);
#endif

//we need to specify explicitly.
    mpreal::set_default_prec(___MPLAPACK_DEFAULT_PRECISION___);
    mpcomplex::set_default_prec(___MPLAPACK_DEFAULT_PRECISION___);

    subst_test1();
    subst_test2();

    addition_real_test();
    subtraction_real_test();
    multiplication_real_test();
    division_real_test();

    addition_complex_test();
    subtraction_complex_test();
    multiplication_complex_test();
    division_complex_test();

    mp_sub_test_real();
    mp_sub_test_complex();

    printf("*** Testing arithmetic successful ***\n");
    return (0);
}
