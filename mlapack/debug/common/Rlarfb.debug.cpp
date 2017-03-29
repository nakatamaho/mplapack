/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlarfb.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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
#include <mblas.h>
#include <mlapack.h>
#include <mpack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_K 1
#define MAX_K 6
#define MIN_N 1
#define MAX_N 6
#define MIN_M 1
#define MAX_M 6
#define MAX_LDT  7
#define MAX_LDC  7
#define MAX_LDV  7
#define MAX_ITER 1

REAL_REF maxdiff = 0.0;

void Rlarfb_test2(const char *side, const char *trans, const char *direct, const char *storev)
{
    int errorflag = FALSE;
    int iter;
    int k, m, n;
    int ldt, ldc;
    int ldv, ldvmin = 0, v_column = 0, ldwork = 0, kmax = 0;
    REAL_REF diff;

    for (n = MIN_N; n <= MAX_N; n++) {
	for (m = MIN_M; m <= MAX_M; m++) {
	    if (Mlsame(side, "L"))
		kmax = m;
	    if (Mlsame(side, "R"))
		kmax = n;
	    for (k = MIN_K; k <= kmax; k++) {

		for (ldt = k; ldt <= MAX_LDT; ldt++) {
		    for (ldc = max(1, m); ldc <= MAX_LDC; ldc++) {
			if (Mlsame(storev, "C") && Mlsame(side, "L"))
			    ldvmin = max(1, m);
			if (Mlsame(storev, "C") && Mlsame(side, "R"))
			    ldvmin = max(1, n);
			if (Mlsame(storev, "R"))
			    ldvmin = k;

			for (ldv = ldvmin; ldv <= MAX_LDV; ldv++) {

			    if (Mlsame(side, "L"))
				ldwork = max(1, n);
			    if (Mlsame(side, "R"))
				ldwork = max(1, m);
			    if (Mlsame(storev, "C"))
				v_column = k;
			    if (Mlsame(storev, "R") && Mlsame(side, "L"))
				v_column = m;
			    if (Mlsame(storev, "R") && Mlsame(side, "R"))
				v_column = n;
#if defined VERBOSE_TEST
			    printf("# side %s: trans %s: direct %s: storev %s, m %d, n %d, k %d, ldv %d, ldt %d, ldc %d, ldwork %d\n",
				 side, trans, direct, storev, m, n, k, ldv, ldt, ldc, ldwork);
#endif
			    REAL_REF *T_ref = new REAL_REF[matlen(ldt, k)];
			    REAL_REF *C_ref = new REAL_REF[matlen(ldc, n)];
			    REAL_REF *V_ref = new REAL_REF[matlen(ldv, v_column)];
			    REAL_REF *work_ref = new REAL_REF[matlen(ldwork, k)];

			    REAL *T = new REAL[matlen(ldt, k)];
			    REAL *C = new REAL[matlen(ldc, n)];
			    REAL *V = new REAL[matlen(ldv, v_column)];
			    REAL *work = new REAL[matlen(ldwork, k)];

			    for (iter = 0; iter < MAX_ITER; iter++) {
				set_random_vector(T_ref, T, matlen(ldt, k));
				set_random_vector(C_ref, C, matlen(ldc, n));
				set_random_vector(V_ref, V, matlen(ldv, v_column));
#if defined ___MPACK_BUILD_WITH_MPFR___
				dlarfb_f77(side, trans, direct, storev, &m, &n, &k, V_ref, &ldv, T_ref, &ldt, C_ref, &ldc, work_ref,
					   &ldwork);
#else
				Rlarfb(side, trans, direct, storev, m, n, k, V_ref, ldv, T_ref, ldt, C_ref, ldc, work_ref, ldwork);
#endif
				Rlarfb(side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc, work, ldwork);

				diff = infnorm(C_ref, C, matlen(ldc, n), 1);
				if (diff > EPSILON) {
				    printf("error: ");  printnum(diff); printf("\n");
				    errorflag = TRUE;
				}
				if (maxdiff < diff)
				    maxdiff = diff;
#if defined VERBOSE_TEST
				printf("max error: "); printnum(maxdiff); printf("\n");
#endif
			    }
			    delete[]work_ref;
			    delete[]V_ref;
			    delete[]C_ref;
			    delete[]T_ref;
			    delete[]work;
			    delete[]V;
			    delete[]C;
			    delete[]T;
			}
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Rlarfb failed ***\n");
	exit(1);
    }
}

void Rlarfb_test()
{
    Rlarfb_test2("L", "N", "F", "C");
    Rlarfb_test2("L", "T", "F", "C");
    Rlarfb_test2("L", "N", "B", "C");
    Rlarfb_test2("L", "T", "B", "C");
    Rlarfb_test2("R", "N", "F", "C");
    Rlarfb_test2("R", "T", "F", "C");
    Rlarfb_test2("R", "N", "B", "C");
    Rlarfb_test2("R", "T", "B", "C");

    Rlarfb_test2("L", "N", "F", "R");
    Rlarfb_test2("L", "T", "F", "R");
    Rlarfb_test2("L", "N", "B", "R");
    Rlarfb_test2("L", "T", "B", "R");
    Rlarfb_test2("R", "N", "F", "R");
    Rlarfb_test2("R", "T", "F", "R");
    Rlarfb_test2("R", "N", "B", "R");
    Rlarfb_test2("R", "T", "B", "R");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rlarfb start ***\n");
    Rlarfb_test();
    printf("*** Testing Rlarfb successful ***\n");
    return (0);
}
