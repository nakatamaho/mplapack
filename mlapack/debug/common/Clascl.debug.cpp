/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clascl.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N     5
#define MAX_N     10
#define MIN_M     5
#define MAX_M     10
#define MIN_LDA   5
#define MAX_LDA   10
#define MAX_ITER   2

REAL_REF maxdiff = 0.0;

void Clascl_test2(const char *type)
{
    int errorflag = FALSE;
    INTEGER_REF j, m, n, lda, ku, kl, info_ref;
    INTEGER info;
    REAL_REF cto_ref, cfrom_ref, diff;
    REAL cto, cfrom;

    for (n = MIN_N; n <= MAX_N; n++) {
	for (m = MIN_M; m <= MAX_M; m++) {
	    for (lda = m; lda <= MAX_LDA; lda++) {
		for (kl = 0; kl <= n; kl++) {
		    for (ku = 0; ku <= m; ku++) {
			COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
			COMPLEX *A = new COMPLEX[matlen(lda, n)];
#if defined VERBOSE_TEST
			printf("#type %s, n %d, m %d, lda %d, ku %d, kl %d\n", type, (int)n, (int)m, (int)lda, (int)ku, (int)kl);
#endif
			j = 0;
			while (j < MAX_ITER) {
			    set_random_vector(A_ref, A, matlen(lda, n));
			    set_random_number(cto_ref, cto);
			    set_random_number(cfrom_ref, cfrom);
#if defined ___MPACK_BUILD_WITH_MPFR___
			    zlascl_f77(type, &kl, &ku, &cfrom_ref, &cto_ref, &m, &n, A_ref, &lda, &info_ref);
#else
			    Clascl(type, kl, ku, cfrom_ref, cto_ref, m, n, A_ref, lda, &info_ref);
#endif
			    Clascl(type, kl, ku, cfrom, cto, m, n, A, lda, &info);

			    if (info_ref != info) {
				printf("error in info %d %d!!\n", (int) info_ref, (int) info);
				errorflag = TRUE;
			    }
                            diff = infnorm(A_ref, A, matlen(lda, n), 1);
			    if (diff > EPSILON10) {
		                printf("error: "); printnum(diff); printf("\n");
		                errorflag = TRUE;
			    }
	                    if (maxdiff < diff)
		                maxdiff = diff;
#if defined VERBOSE_TEST
	                    printf("max error: "); printnum(maxdiff); printf("\n");
#endif
			    j++;
			}
			delete[]A_ref;
			delete[]A;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
        printf("*** Testing Clascl failed ***\n");
	exit(1);
    }
}

void Clascl_test(void)
{
    Clascl_test2("G");
    Clascl_test2("L");
    Clascl_test2("U");
    Clascl_test2("H");
    Clascl_test2("B");
    Clascl_test2("Q");
    Clascl_test2("Z");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clascl start ***\n");
    Clascl_test();
    printf("*** Testing Clascl successful ***\n");
    return (0);
}
