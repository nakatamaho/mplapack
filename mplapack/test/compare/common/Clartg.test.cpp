/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clartg.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

REAL_REF maxdiff = 0.0;

void Clartg_test() {
    int errorflag = FALSE;

    COMPLEX_REF f_ref, g_ref, sn_ref, r_ref;
    REAL_REF cs_ref, diff;
    COMPLEX f, g, sn, r;
    REAL cs;

    int count = 100;

    while (count--) {
        set_random_number(f_ref, f);
        set_random_number(g_ref, g);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlartg_f77(&f_ref, &g_ref, &cs_ref, &sn_ref, &r_ref);
#else
        Clartg(f_ref, g_ref, cs_ref, sn_ref, r_ref);
#endif
        Clartg(f, g, cs, sn, r);

        diff = abs(cs_ref - cs);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        diff = abs(sn_ref - sn);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        diff = abs(r_ref - r);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        f = 0.0;
        f_ref = 0.0;
        set_random_number(g_ref, g);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlartg_f77(&f_ref, &g_ref, &cs_ref, &sn_ref, &r_ref);
#else
        Clartg(f_ref, g_ref, cs_ref, sn_ref, r_ref);
#endif
        Clartg(f, g, cs, sn, r);

        diff = abs(cs_ref - cs);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        diff = abs(sn_ref - sn);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        diff = abs(r_ref - r);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        g = 0.0;
        g_ref = 0.0;
        set_random_number(f_ref, f);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        zlartg_f77(&f_ref, &g_ref, &cs_ref, &sn_ref, &r_ref);
#else
        Clartg(f_ref, g_ref, cs_ref, sn_ref, r_ref);
#endif
        Clartg(f, g, cs, sn, r);

        diff = abs(cs_ref - cs);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        diff = abs(sn_ref - sn);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
        diff = abs(r_ref - r);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif

        if (errorflag == TRUE) {
            printf("*** Testing Clartg failed ***\n");
            exit(1);
        }
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Clartg start ***\n");
    Clartg_test();
    printf("*** Testing Clartg successful ***\n");
    return (0);
}
