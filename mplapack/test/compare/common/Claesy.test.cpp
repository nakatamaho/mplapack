/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Claesy.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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

void Claesy_test() {
    int errorflag = FALSE;
    COMPLEX_REF cs1_ref, rt1_ref, rt2_ref;
    COMPLEX_REF a_ref, b_ref, c_ref, sn1_ref, evscal_ref;
    COMPLEX cs1, rt1, rt2;
    COMPLEX a, b, c, sn1, evscal;

    int count = 100;
    while (count--) {
#if defined VERBOSE_TEST
        printf("Claesy: general random case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);
        set_random_number(c_ref, c);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
        Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
        Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

        if (abs(rt1_ref - rt1) > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (abs(rt2_ref - rt2) > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (abs(cs1_ref - cs1) > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (abs(sn1_ref - sn1) > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (abs(evscal_ref - evscal) > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
/* checking adf = ab case (|a-c| = 2|b|) */
#if defined VERBOSE_TEST
        printf("Claesy: |a-c| = 2|b| case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(c_ref, c);
        b_ref = (a_ref - c_ref) / (REAL_REF)2.0;
        b = (a - c) / (REAL)2.0;

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
        Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
        Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

        if (abs(rt1_ref - rt1) > EPSILON) {
            errorflag = TRUE;
            printf("Error2\n");
        }
        if (abs(rt2_ref - rt2) > EPSILON) {
            errorflag = TRUE;
            printf("Error2\n");
        }
        if (abs(cs1_ref - cs1) > EPSILON) {
            errorflag = TRUE;
            printf("Error2\n");
        }
        if (abs(sn1_ref - sn1) > EPSILON) {
            errorflag = TRUE;
            printf("Error2\n");
        }
        if (abs(evscal_ref - evscal) > EPSILON) {
            errorflag = TRUE;
            printf("Error2\n");
        }
/* checking sm = 0 case */
#if defined VERBOSE_TEST
        printf("Claesy: sm = 0 case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(b_ref, b);
        c_ref = -a_ref;
        c = -a;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
        Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
        Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

        if (abs(rt1_ref - rt1) > EPSILON) {
            errorflag = TRUE;
            printf("Error3\n");
        }
        if (abs(rt2_ref - rt2) > EPSILON) {
            errorflag = TRUE;
            printf("Error3\n");
        }
        if (abs(cs1_ref - cs1) > EPSILON) {
            errorflag = TRUE;
            printf("Error3\n");
        }
        if (abs(sn1_ref - sn1) > EPSILON) {
            errorflag = TRUE;
            printf("Error3\n");
        }
        if (abs(evscal_ref - evscal) > EPSILON) {
            errorflag = TRUE;
            printf("Error3\n");
        }
/*zero eigenvalue case */
#if defined VERBOSE_TEST
        printf("Claesy: zero eigenvalue case\n");
#endif
        set_random_number(a_ref, a);
        b_ref = c_ref = 0.0;
        b = 0.0, c = 0.0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
        Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
        Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

        if (abs(rt1_ref - rt1) > EPSILON) {
            errorflag = TRUE;
            printf("Error4\n");
        }
        if (abs(rt2_ref - rt2) > EPSILON) {
            errorflag = TRUE;
            printf("Error4\n");
        }
        if (abs(cs1_ref - cs1) > EPSILON) {
            errorflag = TRUE;
            printf("Error4\n");
        }
        if (abs(sn1_ref - sn1) > EPSILON) {
            errorflag = TRUE;
            printf("Error4\n");
        }
        if (abs(evscal_ref - evscal) > EPSILON) {
            errorflag = TRUE;
            printf("Error4\n");
        }
/*zero matrix case */
#if defined VERBOSE_TEST
        printf("Claesy: zero matrix case\n");
#endif
        a_ref = b_ref = c_ref = 0.0;
        a = 0.0, b = 0.0, c = 0.0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
        Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
        Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

        if (abs(rt1_ref - rt1) > EPSILON) {
            errorflag = TRUE;
            printf("Error5\n");
        }
        if (abs(rt2_ref - rt2) > EPSILON) {
            errorflag = TRUE;
            printf("Error5\n");
        }
        if (abs(cs1_ref - cs1) > EPSILON) {
            errorflag = TRUE;
            printf("Error5\n");
        }
        if (abs(sn1_ref - sn1) > EPSILON) {
            errorflag = TRUE;
            printf("Error5\n");
        }
        if (abs(evscal_ref - evscal) > EPSILON) {
            errorflag = TRUE;
            printf("Error5\n");
        }
/*Identity matrix case */
#if defined VERBOSE_TEST
        printf("Claesy: identity matrix case\n");
#endif
        set_random_number(a_ref, a);
        set_random_number(c_ref, c);
        b_ref = 0.0;
        b = 0.0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
        Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
        Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

        if (abs(rt1_ref - rt1) > EPSILON) {
            errorflag = TRUE;
            printf("Error6\n");
        }
        if (abs(rt2_ref - rt2) > EPSILON) {
            errorflag = TRUE;
            printf("Error6\n");
        }
        if (abs(cs1_ref - cs1) > EPSILON) {
            errorflag = TRUE;
            printf("Error6\n");
        }
        if (abs(sn1_ref - sn1) > EPSILON) {
            errorflag = TRUE;
            printf("Error6\n");
        }
        if (abs(evscal_ref - evscal) > EPSILON) {
            errorflag = TRUE;
            printf("Error6\n");
        }
        if (errorflag == TRUE) {
            printf("Claesy test failed...\n");
            exit(-1);
        }
    }
/*Identity matrix case */
#if defined VERBOSE_TEST
    printf("Claesy: diagonal matrix case\n");
#endif
    a_ref = c_ref = 1.0;
    b_ref = 0.0;
    a = 1.0, c = 1.0;
    b = 0.0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    zlaesy_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &evscal_ref, &cs1_ref, &sn1_ref);
#else
    Claesy(a_ref, b_ref, c_ref, rt1_ref, rt2_ref, evscal_ref, cs1_ref, sn1_ref);
#endif
    Claesy(a, b, c, rt1, rt2, evscal, cs1, sn1);

    if (abs(rt1_ref - rt1) > EPSILON) {
        errorflag = TRUE;
        printf("Error7\n");
    }
    if (abs(rt2_ref - rt2) > EPSILON) {
        errorflag = TRUE;
        printf("Error7\n");
    }
    if (abs(cs1_ref - cs1) > EPSILON) {
        errorflag = TRUE;
        printf("Error7\n");
    }
    if (abs(sn1_ref - sn1) > EPSILON) {
        errorflag = TRUE;
        printf("Error7\n");
    }
    if (abs(evscal_ref - evscal) > EPSILON) {
        errorflag = TRUE;
        printf("Error8\n");
    }
    if (errorflag == TRUE) {
        printf("*** Testing Claesy failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Claesy start ***\n");
    Claesy_test();
    printf("*** Testing Claesy successful ***\n");
    return (0);
}
