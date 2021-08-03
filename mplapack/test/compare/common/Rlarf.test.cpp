/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlarf.debug.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 2
#define MAX_N 20
#define MIN_M 2
#define MAX_M 20
#define MIN_LDC 2
#define MAX_LDC 20
#define MAX_ITER 2
#define MIN_INCV 1
#define MAX_INCV 3

REAL_REF maxdiff = 0.0;

void Rlarf_test2(const char *side) {
    int errorflag = FALSE;
    INTEGER_REF dimv;
    INTEGER_REF n, m, ldc, incv, iter;
    REAL_REF tau_ref;
    REAL_REF *C_ref, *v_ref, *work_ref;
    REAL_REF diff;
    REAL tau;
    REAL *C, *v, *work;

    for (incv = MIN_INCV; incv <= MAX_INCV; incv++) {
        if (incv == 0)
            continue;
        for (n = MIN_N; n <= MAX_N; n++) {
            for (m = MIN_M; m <= MAX_M; m++) {
                for (ldc = max((INTEGER_REF)1, m); ldc <= MAX_LDC; ldc++) {

                    if (Mlsame(side, "L"))
                        dimv = (1 + (m - 1) * abs(incv));
                    else // side == R
                        dimv = (1 + (n - 1) * abs(incv));

                    C_ref = new REAL_REF[matlen(ldc, n)];
                    v_ref = new REAL_REF[dimv];
                    if (Mlsame(side, "L"))
                        work_ref = new REAL_REF[n];
                    else // side = R
                        work_ref = new REAL_REF[m];

                    C = new REAL[matlen(ldc, n)];
                    v = new REAL[dimv];
                    if (Mlsame(side, "L"))
                        work = new REAL[n];
                    else // side = R
                        work = new REAL[m];

#if defined VERBOSE_TEST
                    printf("# m: %d, n: %d, incv: %d, ldc: %d, side %s\n", (int)m, (int)n, (int)incv, (int)ldc, side);
#endif
                    for (iter = 0; iter < MAX_ITER; iter++) {
                        set_random_vector(C_ref, C, matlen(ldc, n));
                        set_random_vector(v_ref, v, dimv);
                        set_random_number(tau_ref, tau);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                        dlarf_f77(side, &m, &n, v_ref, &incv, &tau_ref, C_ref, &ldc, work_ref);
#else
                        Rlarf(side, m, n, v_ref, incv, tau_ref, C_ref, ldc, work_ref);
#endif
                        Rlarf(side, m, n, v, incv, tau, C, ldc, work);

                        diff = infnorm(C_ref, C, matlen(ldc, n), 1);
                        if (diff > EPSILON) {
                            printf("error: ");
                            printnum(diff);
                            printf("\n");
                            printf("# tau != 0.0\n");
                            printf("# m: %d, n: %d, incv: %d, ldc: %d, side %s\n", (int)m, (int)n, (int)incv, (int)ldc, side);
                            errorflag = TRUE;
                        }
                        if (maxdiff < diff)
                            maxdiff = diff;
#if defined VERBOSE_TEST
                        printf("max error: ");
                        printnum(maxdiff);
                        printf("\n");
#endif
                    }
                    // tau = 0 case
                    for (iter = 0; iter < MAX_ITER; iter++) {
                        set_random_vector(C_ref, C, matlen(ldc, n));
                        set_random_vector(v_ref, v, dimv);
                        tau_ref = 0.0;
                        tau = 0.0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                        dlarf_f77(side, &m, &n, v_ref, &incv, &tau_ref, C_ref, &ldc, work_ref);
#else
                        Rlarf(side, m, n, v_ref, incv, tau_ref, C_ref, ldc, work_ref);
#endif
                        Rlarf(side, m, n, v, incv, tau, C, ldc, work);

                        diff = infnorm(C_ref, C, matlen(ldc, n), 1);
                        if (diff > EPSILON) {
                            printf("error: ");
                            printnum(diff);
                            printf("\n");
                            printf("# tau = 0.0\n");
                            printf("# m: %d, n: %d, incv: %d, ldc: %d, side %s\n", (int)m, (int)n, (int)incv, (int)ldc, side);
                            errorflag = TRUE;
                        }
                        if (maxdiff < diff)
                            maxdiff = diff;
#if defined VERBOSE_TEST
                        printf("max error: ");
                        printnum(maxdiff);
                        printf("\n");
#endif
                    }
                    delete[] work_ref;
                    delete[] v_ref;
                    delete[] C_ref;
                    delete[] work;
                    delete[] v;
                    delete[] C;
                }
            }
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rlarf failed ***\n");
        exit(1);
    }
}

void Rlarf_test() {
    Rlarf_test2("L");
    Rlarf_test2("R");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rlarf start ***\n");
    Rlarf_test();
    printf("*** Testing Rlarf successful ***\n");
    return (0);
}
