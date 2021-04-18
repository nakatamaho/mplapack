/*
 * Copyright (c) 2008-2010
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
#include <mpblas.h>
#include <blas.h>
#include <mplapack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_INCX  -2
#define MAX_INCX   8
#define MIN_INCY  -2
#define MAX_INCY   8
#define MIN_N     -2
#define MAX_N      8
#define MIN_M     -2
#define MAX_M      8
#define MIN_LDA   -2
#define MAX_LDA    8
#define MAX_ITER   2

REAL_REF maxdiff = 0.0;

void Cgeru_test()
{
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	    for (int n = MIN_N; n < MAX_N; n++) {
		for (int m = MIN_M; m < MAX_M; m++) {
		    for (int lda = m; lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
			printf("#m is %d, n is %d, lda is %d, incx is %d, incy is %d\n", m, n, lda, incx, incy);
#endif
			COMPLEX_REF *A_ref;
			COMPLEX_REF *x_ref;
			COMPLEX_REF *y_ref;
			COMPLEX *A;
			COMPLEX *x;
			COMPLEX *y;

			A_ref = new COMPLEX_REF[matlen(lda, n)];
			x_ref = new COMPLEX_REF[veclen(m, incx)];
			y_ref = new COMPLEX_REF[veclen(n, incy)];
			A = new COMPLEX[matlen(lda, n)];
			x = new COMPLEX[veclen(m, incx)];
			y = new COMPLEX[veclen(n, incy)];

			COMPLEX_REF alpha_ref;
			COMPLEX alpha;

			for (int i = 0; i < MAX_ITER; i++) {
			    set_random_number(alpha_ref, alpha);
			    set_random_vector(A_ref, A, matlen(lda, n));
			    set_random_vector(x_ref, x, veclen(m, incx));
			    set_random_vector(y_ref, y, veclen(n, incy));

			    mplapack_errno = 0; blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			    zgeru_f77(&m, &n, &alpha_ref, x_ref, &incx, y_ref, &incy, A_ref, &lda);
			    mplapack_errno1 = blas_errno;
#else
			    Cgeru(m, n, alpha_ref, x_ref, incx, y_ref, incy, A_ref, lda);
			    mplapack_errno1 = mplapack_errno;
#endif
			    Cgeru(m, n, alpha, x, incx, y, incy, A, lda);
			    mplapack_errno2 = mplapack_errno;
#if defined VERBOSE_TEST
			    printf("errno: mplapack %d, ref %d\n", mplapack_errno1, mplapack_errno2);
#endif
			    if (mplapack_errno1 != mplapack_errno2) {
#if defined VERBOSE_TEST
				printf("error in Mxerbla!!\n");
#endif
   	                        errorflag = TRUE;
			    }

			    REAL_REF diff = infnorm(A_ref, A, matlen(lda, n), 1);
			    if (diff > EPSILON) {
#if defined VERBOSE_TEST
				printf("error: "); printnum(diff); printf("\n");
#endif
				errorflag = TRUE;
			    }
			    if (maxdiff < diff)
				maxdiff = diff;
			}
			delete[]y_ref;
			delete[]x_ref;
			delete[]A_ref;
			delete[]y;
			delete[]x;
			delete[]A;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Cgeru failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Cgeru start ***\n");
    Cgeru_test();
    printf("*** Testing Cgeru successful ***\n");
    return (0);
}
