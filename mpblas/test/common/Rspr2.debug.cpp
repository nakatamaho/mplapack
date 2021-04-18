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

#define MIN_INCX   -1
#define MAX_INCX    5
#define MIN_INCY   -1
#define MAX_INCY    5
#define MIN_N      -2
#define MAX_N      10
#define MAX_ITER    3

REAL_REF maxdiff = 0.0;

void Rspr2_test2(const char *uplo)
{
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	    for (int n = MIN_N; n < MAX_N; n++) {
#if defined VERBOSE_TEST
		printf("#n is %d, incx is %d, incy is %d, uplo is %s.\n", n, incx, incy, uplo);
#endif
		REAL_REF alpha_ref;
		REAL_REF *x_ref;
		REAL_REF *y_ref;
		REAL_REF *AP_ref;
		REAL alpha;
		REAL *x;
		REAL *y;
		REAL *AP;

		x_ref = new REAL_REF[veclen(n, incx)];
		y_ref = new REAL_REF[veclen(n, incy)];
		AP_ref = new REAL_REF[vecplen(n)];
		x = new REAL[veclen(n, incx)];
		y = new REAL[veclen(n, incy)];
		AP = new REAL[vecplen(n)];

		for (int i = 0; i < MAX_ITER; i++) {
		    set_random_vector(AP_ref, AP, vecplen(n));
		    set_random_vector(x_ref, x, veclen(n, incx));
		    set_random_vector(y_ref, y, veclen(n, incy));
		    set_random_number(alpha_ref, alpha);

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		    dspr2_f77(uplo, &n, &alpha_ref, x_ref, &incx, y_ref, &incy, AP_ref);
		    mplapack_errno1 = blas_errno;
#else
		    Rspr2(uplo, n, alpha_ref, x_ref, incx, y_ref, incy, AP_ref);
		    mplapack_errno1 = mplapack_errno;
#endif
		    Rspr2(uplo, n, alpha, x, incx, y, incy, AP);
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
		    REAL_REF diff = infnorm(AP_ref, AP, vecplen(n), 1);
		    if (diff > EPSILON) {
#if defined VERBOSE_TEST
			printf("error: "); printnum(diff); printf("\n");
#endif
			errorflag = TRUE;
		    }
		    if (maxdiff < diff)
			maxdiff = diff;
		}
		delete[]AP_ref;
		delete[]y_ref;
		delete[]x_ref;
		delete[]AP;
		delete[]y;
		delete[]x;
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Rspr2 failed ***\n");
	exit(1);
    }
}

void Rspr2_test()
{
    Rspr2_test2("U");
    Rspr2_test2("L");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rspr2 start ***\n");
    Rspr2_test();
    printf("*** Testing Rspr2 successful ***\n");
}
