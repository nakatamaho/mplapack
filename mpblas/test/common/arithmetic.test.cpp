/*
 * Copyright (c) 2008-2026
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
/*
 *
 * This file is a part of MPLAPACK.
 *
 * Basic arithmetic test for REAL/COMPLEX backends.
 *
 */

#include <complex>
#include <mpblas.h>
#include <mplapack_compare_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !defined MAX_ITER
#define MAX_ITER 100
#endif

#if !defined EPSILON
#define EPSILON 1e-14
#endif

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___ || defined ___MPLAPACK_BUILD_WITH_DD___
#include <qd/fpu.h>
static unsigned int g_oldcw;
static void __attribute__((constructor)) mplapack_test_fpu_init(void) { fpu_fix_start(&g_oldcw); }
static void __attribute__((destructor))  mplapack_test_fpu_fini(void) { fpu_fix_end(&g_oldcw); }
#endif

/* Basic deterministic test cases: avoid "edge of overflow/underflow".
 * These values are chosen to exercise signs, zeros, unity, and typical magnitudes.
 */
static const double kRealTestVals[] = {
    0.0,
    1.0,
    -1.0,
    2.5,
    -3.75,
    10.0,
    -10.0,
    0.125,
    3.141592653589793,
    2.718281828459045
};

static const int kNumRealTestVals = (int)(sizeof(kRealTestVals) / sizeof(kRealTestVals[0]));

/* Complex test values are given as (real, imag) pairs. */
static const double kComplexTestVals[][2] = {
    {0.0, 0.0},
    {1.0, 0.0},
    {0.0, 1.0},
    {-1.0, 2.0},
    {2.5, -3.75},
    {3.141592653589793, 2.718281828459045},
    {0.125, 8.0}
};

static const int kNumComplexTestVals = (int)(sizeof(kComplexTestVals) / sizeof(kComplexTestVals[0]));

/* Print-and-abort helpers for deterministic diagnostics. */
static void report_real_failure(const char *label, const REAL_REF &a_ref, const REAL_REF &b_ref, const REAL &a, const REAL &b, const REAL_REF &expected, const REAL &got) {
    REAL_REF diff = abs(expected - got);
    printf("FAILED: %s\n", label);
    printf("a_ref     = ");
    printnum(a_ref);
    printf("\n");
    printf("b_ref     = ");
    printnum(b_ref);
    printf("\n");
    printf("a         = ");
    printnum(a);
    printf("\n");
    printf("b         = ");
    printnum(b);
    printf("\n");
    printf("expected  = ");
    printnum(expected);
    printf("\n");
    printf("got       = ");
    printnum(got);
    printf("\n");
    printf("diff      = ");
    printnum(diff);
    printf("\n");
    exit(1);
}

static void report_complex_failure(const char *label, const COMPLEX_REF &a_ref, const COMPLEX_REF &b_ref, const COMPLEX &a, const COMPLEX &b, const COMPLEX_REF &expected, const COMPLEX &got) {
    REAL_REF diff = abs(expected - got);
    printf("FAILED: %s\n", label);
    printf("a_ref     = ");
    printnum(a_ref);
    printf("\n");
    printf("b_ref     = ");
    printnum(b_ref);
    printf("\n");
    printf("a         = ");
    printnum(a);
    printf("\n");
    printf("b         = ");
    printnum(b);
    printf("\n");
    printf("expected  = ");
    printnum(expected);
    printf("\n");
    printf("got       = ");
    printnum(got);
    printf("\n");
    printf("|diff|    = ");
    printnum(diff);
    printf("\n");
    exit(1);
}

static void check_close_real(const char *label, const REAL_REF &a_ref, const REAL_REF &b_ref, const REAL &a, const REAL &b, const REAL_REF &expected, const REAL &got) {
    REAL_REF diff = abs(expected - got);
    if (abs(diff) > EPSILON) {
        report_real_failure(label, a_ref, b_ref, a, b, expected, got);
    }
}

static void check_close_complex(const char *label, const COMPLEX_REF &a_ref, const COMPLEX_REF &b_ref, const COMPLEX &a, const COMPLEX &b, const COMPLEX_REF &expected, const COMPLEX &got) {
    REAL_REF diff = abs(expected - got);
    if (abs(diff) > EPSILON) {
        report_complex_failure(label, a_ref, b_ref, a, b, expected, got);
    }
}

void subst_test1() {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__], buf3[__MPLAPACK_BUFLEN__];
    memset(buf1, 0, __MPLAPACK_BUFLEN__);
    memset(buf2, 0, __MPLAPACK_BUFLEN__);
    memset(buf3, 0, __MPLAPACK_BUFLEN__);
    REAL tmp1;
    REAL_REF tmp2;

    printf("*** Substitution test 1 ***\n");
    strcpy(buf1, "-1.234567890123456789012345678901234567890123456789012345678901234567890E1");
// tmp1 = buf1;
#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
    #if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
        tmp1 = strtoflt128(buf1, NULL);
    #elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
        tmp1 = strtof128(buf1, NULL);
    #elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
        sscanf(buf1, "%Le", &tmp1);
    #else
        #error "Unsupported MPLAPACK_BINARY128_IO value"
    #endif
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
    sscanf(buf1, "%le", &tmp1);
#elif defined ___MPLAPACK_BUILD_WITH_BINARY80___
    #if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL
        sscanf(buf1, "%Le", &tmp1);
    #elif MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_STRFROMF64X
        tmp1 = strtof64x(buf1, NULL);
    #else
        #error "Unsupported MPLAPACK_BINARY80_IO value"
    #endif
#else
    tmp1 = buf1;
#endif

    tmp2 = tmp1;

    sprintnum(buf2, tmp1);
    sprintnum(buf3, tmp2);

    printf("original  :%s\n", buf1);
    printf("mplib     :%s\n", buf2);
    printf("subst2refm:%s\n", buf3);

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    if (strncmp(buf1, buf2, 19) == 0 && strncmp(buf2, buf3, 19) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
    if (strncmp(buf1, buf2, 65) == 0 && strncmp(buf2, buf3, __MPLAPACK_BUFLEN__) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#elif defined ___MPLAPACK_BUILD_WITH_DD___
    if (strncmp(buf1, buf2, 30) == 0 && strncmp(buf2, buf3, 30) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#elif defined ___MPLAPACK_BUILD_WITH_QD___
    if (strncmp(buf1, buf2, 61) == 0 && strncmp(buf2, buf3, 61) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
    if (strncmp(buf1, buf2, 19) == 0 && strncmp(buf2, buf3, 19) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#elif defined ___MPLAPACK_BUILD_WITH_BINARY128___
    if (strncmp(buf1, buf2, 37) == 0 && strncmp(buf2, buf3, 37) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#elif defined ___MPLAPACK_BUILD_WITH_BINARY80___
    if (strncmp(buf1, buf2, 21) == 0 && strncmp(buf2, buf3, 21) == 0)
        printf("ok!\n");
    else {
        printf("failed!\n");
        exit(1);
    }
#endif
    printf("*** Substitution test 1 successful ***\n");
}

void subst_test2() {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__], buf3[__MPLAPACK_BUFLEN__];
    memset(buf1, 0, __MPLAPACK_BUFLEN__);
    memset(buf2, 0, __MPLAPACK_BUFLEN__);
    memset(buf3, 0, __MPLAPACK_BUFLEN__);

    COMPLEX_REF temp1r, temp2r, temp3r;
    COMPLEX temp1, temp2, temp3;
    REAL_REF diff;

    printf("*** Substitution test2 ***\n");
    temp1r = COMPLEX_REF(0.0, 0.0);
    temp1 = COMPLEX(0.0, 0.0);
    temp2r = COMPLEX_REF(2.0, 0.0);
    temp2 = COMPLEX(2.0, 0.0);

    temp3r = temp1r + temp2r;
    temp3 = temp1 + temp2;

    sprintnum(buf1, temp1r);
    sprintnum(buf2, temp1);
    sprintnum(buf3, temp3);
    printf("temp1r: %s\n", buf1);
    printf("temp1 : %s\n", buf2);
    printf("temp3 : %s\n", buf3);

    printnum(temp1r);
    printf("\n");
    printnum(temp1);
    printf("\n");

    diff = abs(temp1r - temp1);
    printf("diff = ");
    printnum(diff);
    printf("\n");
    if (abs(diff) > EPSILON) {
        printf("*** Substitution test2 failed ***\n");
        exit(1);
    }
    printf("*** Substitution test2 successful ***\n");
}

void mp_rounding2integer() {
    printf("*** Rounding test to integer ***\n");
    REAL_REF temp1r, temp2r;
    REAL temp1, temp2;
    REAL_REF diff;

    set_random_number(temp1r, temp1);
    temp2r = nint(temp1r);
    temp2 = nint(temp1);

    diff = abs(temp2r - temp2);
    printf("diff = ");
    printnum(diff);
    printf("\n");
    if (abs(diff) > EPSILON) {
        printf("*** Rounding test to integer failed ***\n");
        exit(1);
    }
    printf("*** Rounding test to integer successful ***\n");
}

void mp_nint() {
    printf("*** nint test ***\n");
    REAL_REF temp1r, temp2r;
    REAL temp1, temp2;
    REAL_REF diff;

    set_random_number(temp1r, temp1);
    temp2r = nint(temp1r);
    temp2 = nint(temp1);

    diff = abs(temp2r - temp2);
    printf("diff = ");
    printnum(diff);
    printf("\n");
    if (abs(diff) > EPSILON) {
        printf("*** nint test failed ***\n");
        exit(1);
    }
    printf("*** nint test successful ***\n");
}

void mp_sub_test_real() {
    printf("*** REAL sub test ***\n");
    REAL_REF temp1r, temp2r;
    REAL temp1, temp2;
    REAL_REF diff;

    set_random_number(temp1r, temp1);
    set_random_number(temp2r, temp2);

    temp1r = temp1r - temp2r;
    temp1 = temp1 - temp2;

    diff = abs(temp1r - temp1);
    printf("diff = ");
    printnum(diff);
    printf("\n");
    if (abs(diff) > EPSILON) {
        printf("*** REAL sub test failed ***\n");
        exit(1);
    }
    printf("*** REAL sub test successful ***\n");
}

void mp_sub_test_complex() {
    printf("*** COMPLEX sub test ***\n");
    COMPLEX_REF temp1r, temp2r;
    COMPLEX temp1, temp2;
    REAL_REF diff;

    set_random_number(temp1r, temp1);
    set_random_number(temp2r, temp2);

    temp1r = temp1r - temp2r;
    temp1 = temp1 - temp2;

    diff = abs(temp1r - temp1);
    printf("diff = ");
    printnum(diff);
    printf("\n");
    if (abs(diff) > EPSILON) {
        printf("*** COMPLEX sub test failed ***\n");
        exit(1);
    }
    printf("*** COMPLEX sub test successful ***\n");
}

void addition_real_test() {
    printf("*** REAL Addition test ***\n");
    for (int i = 0; i < kNumRealTestVals; ++i) {
        for (int j = 0; j < kNumRealTestVals; ++j) {
            REAL_REF a_ref = kRealTestVals[i];
            REAL_REF b_ref = kRealTestVals[j];
            REAL a = kRealTestVals[i];
            REAL b = kRealTestVals[j];

            // a + b
            REAL_REF expected = a_ref + b_ref;
            REAL got = a + b;
            check_close_real("REAL add: a + b", a_ref, b_ref, a, b, expected, got);

            // a += b
            REAL got2 = a;
            got2 += b;
            check_close_real("REAL add: a += b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        REAL_REF a_ref, b_ref;
        REAL a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        REAL_REF expected = a_ref + b_ref;
        REAL got = a + b;
        check_close_real("REAL add: random a + b", a_ref, b_ref, a, b, expected, got);

        REAL got2 = a;
        got2 += b;
        check_close_real("REAL add: random a += b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** REAL Addition test successful ***\n");
}

void subtraction_real_test() {
    printf("*** REAL Subtraction test ***\n");
    for (int i = 0; i < kNumRealTestVals; ++i) {
        for (int j = 0; j < kNumRealTestVals; ++j) {
            REAL_REF a_ref = kRealTestVals[i];
            REAL_REF b_ref = kRealTestVals[j];
            REAL a = kRealTestVals[i];
            REAL b = kRealTestVals[j];

            // a - b
            REAL_REF expected = a_ref - b_ref;
            REAL got = a - b;
            check_close_real("REAL sub: a - b", a_ref, b_ref, a, b, expected, got);

            // a -= b
            REAL got2 = a;
            got2 -= b;
            check_close_real("REAL sub: a -= b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        REAL_REF a_ref, b_ref;
        REAL a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        REAL_REF expected = a_ref - b_ref;
        REAL got = a - b;
        check_close_real("REAL sub: random a - b", a_ref, b_ref, a, b, expected, got);

        REAL got2 = a;
        got2 -= b;
        check_close_real("REAL sub: random a -= b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** REAL Subtraction test successful ***\n");
}

void multiplication_real_test() {
    printf("*** REAL Multiplication test ***\n");
    for (int i = 0; i < kNumRealTestVals; ++i) {
        for (int j = 0; j < kNumRealTestVals; ++j) {
            REAL_REF a_ref = kRealTestVals[i];
            REAL_REF b_ref = kRealTestVals[j];
            REAL a = kRealTestVals[i];
            REAL b = kRealTestVals[j];

            // a * b
            REAL_REF expected = a_ref * b_ref;
            REAL got = a * b;
            check_close_real("REAL mul: a * b", a_ref, b_ref, a, b, expected, got);

            // a *= b
            REAL got2 = a;
            got2 *= b;
            check_close_real("REAL mul: a *= b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        REAL_REF a_ref, b_ref;
        REAL a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        REAL_REF expected = a_ref * b_ref;
        REAL got = a * b;
        check_close_real("REAL mul: random a * b", a_ref, b_ref, a, b, expected, got);

        REAL got2 = a;
        got2 *= b;
        check_close_real("REAL mul: random a *= b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** REAL Multiplication test successful ***\n");
}

void division_real_test() {
    printf("*** REAL Division test ***\n");
    const REAL_REF denom_min_ref = 0.25;
    const REAL_REF denom_safe_ref = 1.25;
    const REAL denom_min = 0.25;
    const REAL denom_safe = 1.25;    

    for (int i = 0; i < kNumRealTestVals; ++i) {
        for (int j = 0; j < kNumRealTestVals; ++j) {
            REAL_REF a_ref = kRealTestVals[i];
            REAL_REF b_ref = kRealTestVals[j];
            REAL a = kRealTestVals[i];
            REAL b = kRealTestVals[j];

            // Avoid division by 0 (and near-0) in "basic arithmetic" tests.
            if (abs(b_ref) <= denom_min_ref) {
                b_ref = denom_safe_ref;
                b = denom_safe;
            }

            // a / b
            REAL_REF expected = a_ref / b_ref;
            REAL got = a / b;
            check_close_real("REAL div: a / b", a_ref, b_ref, a, b, expected, got);

            // a /= b
            REAL got2 = a;
            got2 /= b;
            check_close_real("REAL div: a /= b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        REAL_REF a_ref, b_ref;
        REAL a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        if (abs(b_ref) <= denom_min_ref) {
            b_ref = denom_safe;
            b = denom_safe;
        }

        REAL_REF expected = a_ref / b_ref;
        REAL got = a / b;
        check_close_real("REAL div: random a / b", a_ref, b_ref, a, b, expected, got);

        REAL got2 = a;
        got2 /= b;
        check_close_real("REAL div: random a /= b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** REAL Division test successful ***\n");
}

void addition_complex_test() {
    printf("*** COMPLEX Addition test ***\n");
    for (int i = 0; i < kNumComplexTestVals; ++i) {
        for (int j = 0; j < kNumComplexTestVals; ++j) {
            COMPLEX_REF a_ref = COMPLEX_REF(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX_REF b_ref = COMPLEX_REF(kComplexTestVals[j][0], kComplexTestVals[j][1]);
            COMPLEX a = COMPLEX(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX b = COMPLEX(kComplexTestVals[j][0], kComplexTestVals[j][1]);

            // a + b
            COMPLEX_REF expected = a_ref + b_ref;
            COMPLEX got = a + b;
            check_close_complex("COMPLEX add: a + b", a_ref, b_ref, a, b, expected, got);

            // a += b
            COMPLEX got2 = a;
            got2 += b;
            check_close_complex("COMPLEX add: a += b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        COMPLEX_REF a_ref, b_ref;
        COMPLEX a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        COMPLEX_REF expected = a_ref + b_ref;
        COMPLEX got = a + b;
        check_close_complex("COMPLEX add: random a + b", a_ref, b_ref, a, b, expected, got);

        COMPLEX got2 = a;
        got2 += b;
        check_close_complex("COMPLEX add: random a += b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** COMPLEX Addition test successful ***\n");
}

void subtraction_complex_test() {
    printf("*** COMPLEX Subtraction test ***\n");
    for (int i = 0; i < kNumComplexTestVals; ++i) {
        for (int j = 0; j < kNumComplexTestVals; ++j) {
            COMPLEX_REF a_ref = COMPLEX_REF(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX_REF b_ref = COMPLEX_REF(kComplexTestVals[j][0], kComplexTestVals[j][1]);
            COMPLEX a = COMPLEX(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX b = COMPLEX(kComplexTestVals[j][0], kComplexTestVals[j][1]);

            // a - b
            COMPLEX_REF expected = a_ref - b_ref;
            COMPLEX got = a - b;
            check_close_complex("COMPLEX sub: a - b", a_ref, b_ref, a, b, expected, got);

            // a -= b
            COMPLEX got2 = a;
            got2 -= b;
            check_close_complex("COMPLEX sub: a -= b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        COMPLEX_REF a_ref, b_ref;
        COMPLEX a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        COMPLEX_REF expected = a_ref - b_ref;
        COMPLEX got = a - b;
        check_close_complex("COMPLEX sub: random a - b", a_ref, b_ref, a, b, expected, got);

        COMPLEX got2 = a;
        got2 -= b;
        check_close_complex("COMPLEX sub: random a -= b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** COMPLEX Subtraction test successful ***\n");
}

void multiplication_complex_test() {
    printf("*** COMPLEX Multiplication test ***\n");
    for (int i = 0; i < kNumComplexTestVals; ++i) {
        for (int j = 0; j < kNumComplexTestVals; ++j) {
            COMPLEX_REF a_ref = COMPLEX_REF(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX_REF b_ref = COMPLEX_REF(kComplexTestVals[j][0], kComplexTestVals[j][1]);
            COMPLEX a = COMPLEX(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX b = COMPLEX(kComplexTestVals[j][0], kComplexTestVals[j][1]);

            // a * b
            COMPLEX_REF expected = a_ref * b_ref;
            COMPLEX got = a * b;
            check_close_complex("COMPLEX mul: a * b", a_ref, b_ref, a, b, expected, got);

            // a *= b
            COMPLEX got2 = a;
            got2 *= b;
            check_close_complex("COMPLEX mul: a *= b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        COMPLEX_REF a_ref, b_ref;
        COMPLEX a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        COMPLEX_REF expected = a_ref * b_ref;
        COMPLEX got = a * b;
        check_close_complex("COMPLEX mul: random a * b", a_ref, b_ref, a, b, expected, got);

        COMPLEX got2 = a;
        got2 *= b;
        check_close_complex("COMPLEX mul: random a *= b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** COMPLEX Multiplication test successful ***\n");
}

void division_complex_test() {
    printf("*** COMPLEX Division test ***\n");
    const REAL_REF denom_min = 0.25;
    const COMPLEX_REF denom_safe_ref = COMPLEX_REF(1.25, 0.75);
    const COMPLEX denom_safe = COMPLEX(1.25, 0.75);

    for (int i = 0; i < kNumComplexTestVals; ++i) {
        for (int j = 0; j < kNumComplexTestVals; ++j) {
            COMPLEX_REF a_ref = COMPLEX_REF(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX_REF b_ref = COMPLEX_REF(kComplexTestVals[j][0], kComplexTestVals[j][1]);
            COMPLEX a = COMPLEX(kComplexTestVals[i][0], kComplexTestVals[i][1]);
            COMPLEX b = COMPLEX(kComplexTestVals[j][0], kComplexTestVals[j][1]);

            // Avoid division by 0 (and near-0) in "basic arithmetic" tests.
            if (abs(b_ref) <= denom_min) {
                b_ref = denom_safe_ref;
                b = denom_safe;
            }

            // a / b
            COMPLEX_REF expected = a_ref / b_ref;
            COMPLEX got = a / b;
            check_close_complex("COMPLEX div: a / b", a_ref, b_ref, a, b, expected, got);

            // a /= b
            COMPLEX got2 = a;
            got2 /= b;
            check_close_complex("COMPLEX div: a /= b", a_ref, b_ref, a, b, expected, got2);
        }
    }

    // Random smoke tests
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        COMPLEX_REF a_ref, b_ref;
        COMPLEX a, b;
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);

        if (abs(b_ref) <= denom_min) {
            b_ref = denom_safe_ref;
            b = denom_safe;
        }

        COMPLEX_REF expected = a_ref / b_ref;
        COMPLEX got = a / b;
        check_close_complex("COMPLEX div: random a / b", a_ref, b_ref, a, b, expected, got);

        COMPLEX got2 = a;
        got2 /= b;
        check_close_complex("COMPLEX div: random a /= b", a_ref, b_ref, a, b, expected, got2);
    }

    printf("*** COMPLEX Division test successful ***\n");
}

int main() {
    printf("*** Testing arithmetic ***\n");

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    mpf_set_default_prec(___MPLAPACK_GMP_DEFAULT_PRECISION___);
#endif

    // we need to specify explicitly.
    mpfr_class::default_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpc_class::default_real_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpc_class::default_imag_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;

    mp_rounding2integer();
    mp_nint();

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
