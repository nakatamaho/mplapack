/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlamch.debug.cpp,v 1.11 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Rlamch.debug.cpp,v 1.11 2010/08/07 05:50:10 nakatamaho Exp $

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

#include <mpblas.h>
#include <mplapack.h>
#include <blas.h>
#include <lapack.h>
#include <mplapack_debug.h>
#include <iostream>

#define VERBOSE_TEST

/* dlamch result on FreeBSD 6/i386 + Lapack 3.1.1
  Epsilon                      =   1.110223024625157E-016
  Safe minimum                 =   2.225073858507201E-308
  Base                         =    2.00000000000000     
  Precision                    =   2.220446049250313E-016
  Number of digits in mantissa =    53.0000000000000     
  Rounding mode                =    1.00000000000000     
  Minimum exponent             =   -1021.00000000000     
  Underflow threshold          =   2.225073858507201E-308
  Largest exponent             =    1024.00000000000     
  Overflow threshold           =   1.797693134862316E+308
*/

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
void Rlamch_mpfr_test()
{
    REAL tmp, tmp2;
    tmp = Rlamch_mpfr("E");
    tmp2 = Rlamch_mpfr("E");
#if defined VERBOSE_TEST
    printf("Rlamch E: Epsilon                      "); printnum(Rlamch_mpfr("E")); printf("\n");
#endif

    if (1.0 + tmp * 2.0 <= 1.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }

    tmp = Rlamch_mpfr("S");
    printf("Rlamch S: Safe minimum                 "); printnum(tmp); printf("\n");

/*
    REAL tmp3 = tmp / 2.0;
    REAL tmp4 = 1.0 / (tmp);
    REAL tmp5 = 1.0 / (tmp3);
    REAL tmp6 = 2.225073858507201e-308;
    mpfr_printf("%72.64Rb\n", mpfr_ptr(tmp));
    mpfr_printf("%72.64Rb should be smaller than above\n", mpfr_ptr(tmp3)); 
    mpfr_printf("%72.64Rb\n", mpfr_ptr(tmp4));
    mpfr_printf("%72.64Rb should overflow\n", mpfr_ptr(tmp5));
    mpfr_printf("%72.64Rb just a ref.\n", mpfr_ptr(tmp6));
*/

    tmp = Rlamch_mpfr("B");
    printf("Rlamch B: Base                         "); printnum(tmp); printf("\n");

    tmp = Rlamch_mpfr("P");
    printf("Rlamch P: Precision                    "); printnum(tmp); printf("\n");

    tmp = Rlamch_mpfr("N");
    printf("Rlamch N: Number of digits in mantissa "); printnum(tmp); printf("\n");

    tmp = Rlamch_mpfr("R");
    printf("Rlamch R: Rounding mode                "); printnum(tmp); printf("\n");

    tmp = Rlamch_mpfr("M");
    printf("Rlamch M: Minimum exponent:            "); printnum(tmp); printf("\n");

/*
    REAL tmp3, tmp4, tmp5;
    REAL one = 1.0;
    tmp3= mul_2si(one, (long)tmp);
    tmp4= (1.0 - Rlamch_mpfr("E")) * mul_2si(one, (long)tmp); //fill out sig. dig. by 1.
    tmp5= tmp4 / 2.0;
    tmp5= tmp5 * 2.0;
    mpfr_printf("%72.64Rb 2^(minexp)\n", mpfr_ptr(tmp3));
    mpfr_printf("%72.64Rb fill out by 1\n", mpfr_ptr(tmp4));
    mpfr_printf("%72.64Rb abrupt undeflow MPFR; gradual underflow: lost one bit\n", mpfr_ptr(tmp5));
*/
//   
    tmp = Rlamch_mpfr("U");
    printf("Rlamch U: Underflow threshold          "); printnum(tmp); printf("\n");
/*
    mpfr_printf("%72.64Rb\n", mpfr_ptr(tmp));
*/

    tmp = Rlamch_mpfr("L");
    printf("Rlamch L: Largest exponent             "); printnum(tmp); printf("\n");

/*
    REAL tmp3, tmp4;
    REAL one = 1.0;
    tmp3= mul_2si(one, (long)tmp);
    tmp4 = tmp3*2.0;
    mpfr_printf("%72.64Rb 2^(largeexp)\n", mpfr_ptr(tmp3));
    mpfr_printf("%72.64Rb should overflow\n", mpfr_ptr(tmp4));
*/

    tmp = Rlamch_mpfr("O");
    printf("Rlamch O: Overflow threshold           "); printnum(tmp); printf("\n");
/*
    mpfr_printf("%72.64Rb over thr.\n", mpfr_ptr(tmp));
*/
    tmp = Rlamch_mpfr("S");
    tmp = 1.0 / tmp;
    printf("Rlamch -: Reciprocal of safe minimum   "); printnum(tmp); printf("\n");

    tmp = Rlamch_mpfr("Z");
    printf("Rlamch Z: dummy (error)                "); printnum(tmp); printf("\n");
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void Rlamch_gmp_test()
{
    REAL tmp, tmp2;
    tmp = Rlamch_gmp("E");
    tmp2 = Rlamch_gmp("E");
#if defined VERBOSE_TEST
    printf("Rlamch E: Epsilon                      "); printnum(Rlamch_gmp("E")); printf("\n");
#endif
    if (abs(tmp - tmp2) != 0.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }
    if (1.0 + tmp*0.5 != 1.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }
    
    tmp = Rlamch_gmp("S");
    printf("Rlamch S: Safe minimum                 "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("B");
    printf("Rlamch B: Base                         "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("P");
    printf("Rlamch P: Precision                    "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("N");
    printf("Rlamch N: Number of digits in mantissa "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("R");
    printf("Rlamch R: Rounding mode                "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("M");
    printf("Rlamch M: Minimum exponent:            "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("U");
    printf("Rlamch U: Underflow threshold          "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("L");
    printf("Rlamch L: Largest exponent             "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("O");
    printf("Rlamch O: Overflow threshold           "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("S");
    tmp = 1.0 / tmp;
    printf("Rlamch -: Reciprocal of safe minimum   "); printnum(tmp); printf("\n");

    tmp = Rlamch_gmp("Z");
    printf("Rlamch Z: dummy (error)                "); printnum(tmp); printf("\n");
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
void
Rlamch_qd_test()
{
    qd_real tmp, tmp2;
    tmp = Rlamch_qd("E");
    tmp2 = Rlamch_qd("E");
#if defined VERBOSE_TEST
    printf("Rlamch E: Epsilon                      "); printnum(Rlamch_qd("E")); printf("\n");
#endif
    if (abs(tmp - tmp2) != 0.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }
    if (1.0 + tmp <= 1.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }
    tmp = Rlamch_qd("S");
    printf("Rlamch S: Safe minimum                 "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("B");
    printf("Rlamch B: Base                         "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("P");
    printf("Rlamch P: Precision                    "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("N");
    printf("Rlamch N: Number of digits in mantissa "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("R");
    printf("Rlamch R: Rounding mode                "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("M");
    printf("Rlamch M: Minimum exponent:            "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("U");
    printf("Rlamch U: Underflow threshold          "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("L");
    printf("Rlamch L: Largest exponent             "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("O");
    printf("Rlamch O: Overflow threshold           "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("S");
    tmp = 1.0 / tmp;
    printf("Rlamch -: Reciprocal of safe minimum   "); printnum(tmp); printf("\n");

    tmp = Rlamch_qd("Z");
    printf("Rlamch Z: dummy (error)                "); printnum(tmp); printf("\n");
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
void
Rlamch_dd_test()
{
    dd_real tmp, tmp2;
    tmp = Rlamch_dd("E");
    tmp2 = Rlamch_dd("E");
#if defined VERBOSE_TEST
    printf("Rlamch E: Epsilon                      "); printnum(Rlamch_dd("E")); printf("\n");
#endif
    if (abs(tmp - tmp2) != 0.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }
    if (1.0 + tmp <= 1.0) {
        printf("*** Testing Mutils (RlamchE) failed ***\n");
	exit(1);
    }
    tmp = Rlamch_dd("S");
    printf("Rlamch S: Safe minimum                 "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("B");
    printf("Rlamch B: Base                         "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("P");
    printf("Rlamch P: Precision                    "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("N");
    printf("Rlamch N: Number of digits in mantissa "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("R");
    printf("Rlamch R: Rounding mode                "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("M");
    printf("Rlamch M: Minimum exponent:            "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("U");
    printf("Rlamch U: Underflow threshold          "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("L");
    printf("Rlamch L: Largest exponent             "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("O");
    printf("Rlamch O: Overflow threshold           "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("S");
    tmp = 1.0 / tmp;
    printf("Rlamch -: Reciprocal of safe minimum   "); printnum(tmp); printf("\n");

    tmp = Rlamch_dd("Z");
    printf("Rlamch Z: dummy (error)                "); printnum(tmp); printf("\n");
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
void Rlamch_double_test()
{
    double tmp;
    tmp = Rlamch_double("E");
#if defined VERBOSE_TEST
    printf("Rlamch E: Epsilon                      "); printnum(Rlamch_double("E")); printf("\n");
#endif
    tmp = Rlamch_double("S");
    printf("Rlamch S: Safe minimum                 "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("B");
    printf("Rlamch B: Base                         "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("P");
    printf("Rlamch P: Precision                    "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("N");
    printf("Rlamch N: Number of digits in mantissa "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("R");
    printf("Rlamch R: Rounding mode                "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("M");
    printf("Rlamch M: Minimum exponent:            "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("U");
    printf("Rlamch U: Underflow threshold          "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("L");
    printf("Rlamch L: Largest exponent             "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("O");
    printf("Rlamch O: Overflow threshold           "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("S");
    tmp = 1.0 / tmp;
    printf("Rlamch -: Reciprocal of safe minimum   "); printnum(tmp); printf("\n");

    tmp = Rlamch_double("Z");
    printf("Rlamch Z: dummy (error)                "); printnum(tmp); printf("\n");
}
#endif

#if defined ___MPLAPACK_BUILD_WITH___FLOAT128___
void Rlamch___float128_test()
{
    __float128 tmp;
    tmp = Rlamch___float128("E");
/*
    if (1.0Q + tmp > 1.0Q) {printf("rlamche f128 ok\n");} else {printf("rlamche f128 Error\n");}
    if (1.0Q + tmp/2.0Q > 1.0Q) {printf("rlamche f128 error\n");} else {printf("rlamche f128 ok\n");}
    if (1.0Q + tmp/2.0Q == 1.0Q) {printf("rlamche f128 ok\n");} else {printf("rlamche f128 error\n");}
*/
#if defined VERBOSE_TEST
    printf("Rlamch E: Epsilon                      "); printnum(Rlamch___float128("E")); printf("\n");
#endif
    tmp = Rlamch___float128("S");
    printf("Rlamch S: Safe minimum                 "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("B");
    printf("Rlamch B: Base                         "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("P");
    printf("Rlamch P: Precision                    "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("N");
    printf("Rlamch N: Number of digits in mantissa "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("R");
    printf("Rlamch R: Rounding mode                "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("M");
    printf("Rlamch M: Minimum exponent:            "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("U");
    printf("Rlamch U: Underflow threshold          "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("L");
    printf("Rlamch L: Largest exponent             "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("O");
    printf("Rlamch O: Overflow threshold           "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("S");
    tmp = 1.0 / tmp;
    printf("Rlamch -: Reciprocal of safe minimum   "); printnum(tmp); printf("\n");

    tmp = Rlamch___float128("Z");
    printf("Rlamch Z: dummy (error)                "); printnum(tmp); printf("\n");
}
#endif

int main(int argc, char *argv[])
{
    printf("*** Testing Rlamch start ***\n");
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    Rlamch_mpfr_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_GMP___
    Rlamch_gmp_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_QD___
    Rlamch_qd_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_DD___
    Rlamch_dd_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
    Rlamch_double_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH___FLOAT128___
    Rlamch___float128_test();
#endif
    printf("*** Testing Rlamch successful ***\n");
    return(0);
}
