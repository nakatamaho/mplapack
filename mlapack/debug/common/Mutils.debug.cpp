/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Mutils.debug.cpp,v 1.13 2010/08/07 05:50:10 nakatamaho Exp $
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *
 * $Id: Mutils.debug.cpp,v 1.13 2010/08/07 05:50:10 nakatamaho Exp $

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mblas.h>
#include <mlapack.h>
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#include <cmath>

using namespace std;

#define ITERATION 100
#define MAXSIZE 10
#define MAXCOND 16

void Mutils_test_pi()
{
    int errorflag = TRUE;
    REAL_REF p_ref, diff, dummy_ref = 0.0;
    REAL p, dummy = 0.0;

    p_ref = pi(dummy_ref);
    p = pi(dummy);
#if defined VERBOSE_TEST
    cout << "p=       "; printnum(p); cout << endl;
    cout << "p_ref=   "; printnum(p_ref); cout << endl;
    cout << "residue=p-p_ref" << endl;
#endif
    diff = abs( p_ref - p);
#if defined VERBOSE_TEST
    printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	errorflag = TRUE;
        printf("*** Testing Mutils (pi) failed ***\n");
	exit(1);
    }
}

//ln2 test
void Mutils_test_log2()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	a_ref = abs(a_ref);
        a = abs(a);

	b_ref = log2(a_ref);
#if defined ___MPACK_BUILD_WITH___FLOAT128___
        b = std::log2(a);
#else
        b = log2(a);
#endif


#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "ln2_ref= "; printnum(b_ref); cout << endl;
	cout << "ln2a=    "; printnum(b); cout << endl;
#endif
        diff = abs( b_ref - b);
#if defined VERBOSE_TEST
        printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (log2) failed ***\n");
	    exit(1);
	}
    }
}

//ln test
void Mutils_test_log()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	a_ref = abs(a_ref);
	a = abs(a);

	b_ref = log(a_ref);
	b = log(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "lna_ref= "; printnum(b_ref); cout << endl;
	cout << "lna=     "; printnum(b); cout << endl;
#endif
        diff = abs (b_ref - b);
#if defined VERBOSE_TEST
        printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (log) failed ***\n");
	    exit(1);
	}
    }
}

void Mutils_test_log10()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	a_ref = abs(a_ref);
	b_ref = log10(a_ref);
        a = abs(a);
	b = log10(a);
#if defined VERBOSE_TEST
	cout << "a_ref=      "; printnum(a_ref); cout << endl;
	cout << "a=          "; printnum(a); cout << endl;
	cout << "log10a_ref= "; printnum(b_ref);cout << endl;
	cout << "log10a=     "; printnum(b); cout << endl;
#endif
        diff = abs (b_ref - b);
#if defined VERBOSE_TEST
        printf("diff        "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (log10) failed ***\n");
	    exit(1);
	}
    }
}

void Mutils_test_sign()
{
    int errorflag = TRUE;
    REAL mtemp;

    mtemp = sign((REAL)1.0, (REAL)-1.0);
    if (mtemp != -1.0) {
        printf("*** Testing Mutils (sign) failed ***\n");
	errorflag = FALSE;
        exit(1);
    }
    mtemp = sign((REAL)-1.0, (REAL)1.0);
    if (mtemp != 1.0) {
        printf("*** Testing Mutils (sign) failed ***\n");
	errorflag = FALSE;
        exit(1);
    }
    mtemp = sign((REAL)-1.0, (REAL)0.0);
    if (mtemp != 1.0) {
        printf("*** Testing Mutils (sign) failed ***\n");
	errorflag = FALSE;
        exit(1);
    }
    mtemp = sign((REAL)1.0, (REAL)0.0);
    if (mtemp != 1.0) {
        printf("*** Testing Mutils (sign) failed ***\n");
	errorflag = FALSE;
        exit(1);
    }

    mtemp = sign((REAL)0.0, (REAL)0.0);
    if (mtemp != 0.0) {
        printf("*** Testing Mutils (sign) failed ***\n");
	errorflag = FALSE;
        exit(1);
    }
    if (errorflag == FALSE) {
        printf("*** Testing Mutils (sign) failed ***\n");
	errorflag = FALSE;
        exit(1);
    }
}

void Mutils_test_pow()
{
    int errorflag = TRUE;
    REAL_REF x_ref, y_ref, z_ref, diff;
    REAL x, y, z;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(x_ref, x);
	set_random_number(y_ref, y);
	x_ref = abs(x_ref);
	y_ref = abs(y_ref);
#if defined ___MPACK_BUILD_WITH_MPFR___
        z_ref = std::pow(x_ref, y_ref);
#else
        z_ref = pow(x_ref, y_ref);
#endif
	x = abs(x);
	y = abs(y);
#if defined ___MPACK_BUILD_WITH_DOUBLE___
        z = std::pow(x, y);
#else
	z = pow(x, y);
#endif

#if defined VERBOSE_TEST
	cout << "x_ref=   "; printnum(x_ref); cout << endl;
	cout << "y_ref=   "; printnum(y_ref); cout << endl;
	cout << "z_ref=   "; printnum(z_ref); cout << endl;
	cout << "x=       "; printnum(x); cout << endl;
	cout << "y=       "; printnum(y); cout << endl;
	cout << "z=       "; printnum(z); cout << endl;
#endif
        diff = abs (z_ref - z);
#if defined VERBOSE_TEST
        printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (pow) failed ***\n");
	    exit(1);
	}
    }
}

//sin test
void Mutils_test_sin()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	b_ref = sin(a_ref);
	b = sin(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "sin_ref= "; printnum(b_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "sin=     "; printnum(b); cout << endl;
	cout << endl;
#endif
        diff = abs( b_ref - b);
#if defined VERBOSE_TEST
        printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (sin) failed ***\n");
	    exit(1);
	}
    }
}

//cos test
void Mutils_test_cos()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	b_ref = cos(a_ref);
	b = cos(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "cos_ref= "; printnum(b_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "cos=     "; printnum(b); cout << endl;
#endif
        diff = abs( b_ref - b);
#if defined VERBOSE_TEST
        printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (cos) failed ***\n");
	    exit(1);
	}
    }
}

//exp test
void Mutils_test_exp()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	b_ref = exp(a_ref);
	b = exp(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "exp_ref= "; printnum(b_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "exp=     "; printnum(b); cout << endl;
#endif
        diff = abs( b_ref - b);
#if defined VERBOSE_TEST
        printf("diff     "); printnum(diff); printf("\n\n");
#endif
#if defined ___MPACK_BUILD_WITH_GMP___
    if (diff > EPSILON100) {
#else
    if (diff > EPSILON) {
#endif
	    errorflag = TRUE;
            printf("*** Testing Mutils (exp) failed ***\n");
	    exit(1);
	}
    }
}

void Mutils_test_highcond()
{
  for (int n = 2 ; n < MAXSIZE;  n++ ) {
    REAL_REF *A_ref = new REAL_REF [n * n];
    REAL *A = new REAL [n * n];

    for (int c = 1 ; c < MAXCOND;  c++ ) {
      set_random_symmmat_cond(A_ref, A, n, n, c); 
      printf("approx cond: %d\n",c);
      printf("A="); printmat(n, n, A, n); printf("\n");printf("\n");
    }
    delete []A;
    delete []A_ref;
  }
}

void Mutils_test()
{
    Mutils_test_pi();
    Mutils_test_log2();
    Mutils_test_log();
    Mutils_test_log10();
    Mutils_test_sign();
    Mutils_test_pow();
    Mutils_test_sin();
    Mutils_test_cos();
    Mutils_test_exp();
//  Mutils_test_highcond();
}

int main(int argc, char *argv[])
{
    printf("*** Testing Mutils start ***\n");
    Mutils_test();
    printf("*** Testing Mutils successful ***\n");
    return (0);
}
