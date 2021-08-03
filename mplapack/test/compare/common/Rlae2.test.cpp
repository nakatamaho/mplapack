/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlae2.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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

#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_compare_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define ITERATION 100

REAL_REF maxdiff = 0.0;

void Rlae2_test() {
    int errorflag = FALSE;
    REAL_REF a_ref, b_ref, c_ref, rt1_ref, rt2_ref, diff;
    REAL a, b, c, rt1, rt2;
    int count = ITERATION;

    while (count--) {
#if defined VERBOSE_TEST
        printf("Rlae2: general random case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);
        set_random_number(c_ref, c);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        dlae2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref);
#else
        Rlae2(a_ref, b_ref, c_ref, rt1_ref, rt2_ref);
#endif
        Rlae2(a, b, c, rt1, rt2);

        diff = abs(rt1_ref - rt1);
#if defined VERBOSE_TEST
        printf("diff1    ");
        printnum(diff);
        printf("\n");
#endif
        if (diff > EPSILON) {
            errorflag = TRUE;
        }
        if (maxdiff < diff)
            maxdiff = diff;
        diff = abs(rt2_ref - rt2);
#if defined VERBOSE_TEST
        printf("diff1    ");
        printnum(diff);
        printf("\n");
#endif
        if (diff > EPSILON) {
            errorflag = TRUE;
        }
        if (maxdiff < diff)
            maxdiff = diff;

/* checking adf = ab (= O) case (|a-c| = 2|b|) */
#if defined VERBOSE_TEST
        printf("Rlae2: |a-c|=2|b| case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(c_ref, c);
        b_ref = (a_ref - c_ref) * 0.5;
        b = (a - c) * 0.5;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        dlae2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref);
#else
        Rlae2(a_ref, b_ref, c_ref, rt1_ref, rt2_ref);
#endif
        Rlae2(a, b, c, rt1, rt2);

        diff = abs(rt1_ref - rt1);
#if defined VERBOSE_TEST
        printf("diff1    ");
        printnum(diff);
        printf("\n");
#endif
        if (diff > EPSILON) {
            errorflag = TRUE;
        }
        if (maxdiff < diff)
            maxdiff = diff;

        diff = abs(rt2_ref - rt2);
#if defined VERBOSE_TEST
        printf("diff1    ");
        printnum(diff);
        printf("\n");
#endif
        if (diff > EPSILON) {
            errorflag = TRUE;
        }
        if (maxdiff < diff)
            maxdiff = diff;

/* checking rt1=rt2 = 0 case */
#if defined VERBOSE_TEST
        printf("Rlae2: rt1=rt2=0 case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);
        c_ref = -a_ref;
        c = -a;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        dlae2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref);
#else
        Rlae2(a_ref, b_ref, c_ref, rt1_ref, rt2_ref);
#endif
        Rlae2(a, b, c, rt1, rt2);

        diff = abs(rt1_ref - rt1);
#if defined VERBOSE_TEST
        printf("diff1    ");
        printnum(diff);
        printf("\n");
#endif
        if (diff > EPSILON) {
            errorflag = TRUE;
        }
        if (maxdiff < diff)
            maxdiff = diff;

        diff = abs(rt2_ref - rt2);
#if defined VERBOSE_TEST
        printf("diff1    ");
        printnum(diff);
        printf("\n");
#endif
        if (diff > EPSILON) {
            errorflag = TRUE;
        }
        if (maxdiff < diff)
            maxdiff = diff;

        if (errorflag == TRUE) {
            printf("*** Testing Rlae2 failed ***\n");
            exit(1);
        }
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rlae2 start ***\n");
    Rlae2_test();
    printf("*** Testing Rlae2 successful ***\n");
    return (0);
}
