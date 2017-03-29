/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cungtr.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 0
#define MAX_N 15
#define MIN_M 0
#define MAX_M 15
#define MIN_K 0
#define MAX_K 15
#define MIN_LDA 0
#define MAX_LDA 15
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Cungtr_test2(const char *uplo)
{
    int errorflag = FALSE;
    int iter;
    int n, lda, lwork;
    REAL_REF diff;
    INTEGER_REF info_ref, worksize_ref;
    INTEGER info, worksize;

    for (n = MIN_N; n <= MAX_N; n++) {
	for (lda = max(1, n); lda <= MAX_LDA; lda++) {
#if defined VERBOSE_TEST
	    printf("# uplo %s, n %d, lda %d\n", uplo, n, lda);
#endif
	    COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
	    COMPLEX_REF *tau_ref = new COMPLEX_REF[veclen(n - 1, 1)];
	    COMPLEX_REF *work_ref = new COMPLEX_REF[veclen(n - 1, 1) * 1024];

	    COMPLEX *A = new COMPLEX[matlen(lda, n)];
	    COMPLEX *tau = new COMPLEX[veclen(n - 1, 1)];
	    COMPLEX *work = new COMPLEX[veclen(n - 1, 1) * 1024];

//these workspace query might not be the same value.
	    lwork = -1;
#if defined ___MPACK_BUILD_WITH_MPFR___
	    zungtr_f77(uplo, &n, A_ref, &lda, tau_ref, work_ref, &lwork, &info_ref);
#else
	    Cungtr(uplo, n, A_ref, lda, tau_ref, work_ref, lwork, &info_ref);
#endif
	    Cungtr(uplo, n, A, lda, tau, work, lwork, &info);

	    worksize_ref = (int) cast2double(work_ref[0].real());
	    worksize = (int) cast2double(work[0].real());

#if defined VERBOSE_TEST
	    printf("optimized worksize by dorgtr %d : by Cungtr %d.\n", (int)worksize_ref, (int)worksize);
#endif
#ifdef DUMMY
//comparison of workspace is nonsense...
	    if (worksize != worksized)
		printf("error in worksize\n");
#endif
	    for (iter = 0; iter < MAX_ITER; iter++) {
		set_random_vector(A_ref, A, matlen(lda, n));
		set_random_vector(tau_ref, tau, veclen(n - 1, 1));
		set_random_vector(work_ref, work, veclen(n - 1, 1) * 1024);

		lwork = worksize_ref;
#if defined ___MPACK_BUILD_WITH_MPFR___
		zungtr_f77(uplo, &n, A_ref, &lda, tau_ref, work_ref, &lwork, &info_ref);
#else
		Cungtr(uplo, n, A_ref, lda, tau_ref, work_ref, lwork, &info_ref);
#endif
		Cungtr(uplo, n, A, lda, tau, work, lwork, &info);

		diff = infnorm(A_ref, A, matlen(lda, n), 1);
		if (diff > EPSILON2) {
		    printf("error in A: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff)
		    maxdiff = diff;

		diff = infnorm(tau_ref, tau, veclen(n - 1, 1), 1);
		if (diff > EPSILON2) {
		    printf("error in tau: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff)
		    maxdiff = diff;
#ifdef DUMMY
//comparison of workspace is nonsense...
		diff = infnorm(work_ref, work, veclen(n - 1, 1), 1);
		if (diff > EPSILON) {
		    printf("error in work: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff)
		    maxdiff = diff;
#endif
#if defined VERBOSE_TEST
   	        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
	    }
	    delete[]tau_ref;
	    delete[]work_ref;
	    delete[]A_ref;
	    delete[]tau;
	    delete[]work;
	    delete[]A;
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Cungtr failed ***\n");
	exit(1);
    }
}

void Cungtr_test()
{
    Cungtr_test2("L");
    Cungtr_test2("U");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Cungtr start ***\n");
    Cungtr_test();
    printf("*** Testing Cungtr successful ***\n");
    return (0);
}
