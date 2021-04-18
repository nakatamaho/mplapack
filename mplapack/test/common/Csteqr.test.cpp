/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Csteqr.debug.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
#include <mplapack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N      2
#define MAX_N     20
#define MAX_LDZ   20
#define MAX_ITER   2

template < class X, class Y > X infnorm(X * vec1, Y * vec2, int len, int inc)
{
    Y ctemp;
    X inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
	ctemp = vec1[i] - vec2[i];
	if (inorm < abs(ctemp)) {
	    inorm = abs(ctemp);
	}
    }
    return inorm;
}

template < class X, class Y > X check_eigvec(X * d, X * eig, X * e, Y * Z, int n, int ldz)
{
    Y *tmpmat1 = new Y[n * n];	//original mat 
    X *tmpmat2 = new X[n * n];	//eigen mat
    Y *tmpmat3 = new Y[n * n];
    X diff;
    Y mtemp;

    for (int i = 0; i < n * n; i++) {
	tmpmat1[i] = 0.0;
	tmpmat2[i] = 0.0;
    }

    for (int i = 0; i < n; i++) {
	tmpmat1[i + i * n] = d[i];
	tmpmat2[i + i * n] = eig[i];
    }
    for (int i = 1; i < n; i++) {
	tmpmat1[i + (i - 1) * n] = e[i - 1];
	tmpmat1[(i - 1) + i * n] = e[i - 1];
    }
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtemp = 0.0;
	    for (int p = 0; p < n; p++) {
		for (int q = 0; q < n; q++) {
		    mtemp = mtemp + Z[p + i * ldz] * tmpmat1[p + q * n] * Z[q + j * ldz];
		}
	    }
	    tmpmat3[i + j * n] = mtemp;
	}
    }
    diff = infnorm < X, Y > (tmpmat2, tmpmat3, n * n, 1);
    return diff;
/*
  ( Z' * X * Z ) = Z_
  printf("X = ");  printmat(n, n, tmpmat1, n); printf("\n\n");
  printf("E = ");  printmat(n, n, tmpmat2, n); printf("\n\n");
  printf("Z_= ");  printmat(n, n, tmpmat3, n); printf("\n\n");
  printf("Z = ");  printmat(n, n, Z, ldz); printf("\n\n");
*/
}

REAL_REF maxdiff = 0.0;

void Csteqr_test2(const char *compz)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, ldz_ref, n_ref;
    INTEGER info;
    REAL_REF diff = 0.0;

    for (INTEGER n = MIN_N; n < MAX_N; n++) {
	for (INTEGER ldz = max(n, (INTEGER) 1); ldz < MAX_LDZ; ldz++) {
	    COMPLEX_REF *Z_ref = new COMPLEX_REF[matlen(ldz, n)];
	    REAL_REF *D_ref = new REAL_REF[veclen(n, 1)];
	    REAL_REF *E_ref = new REAL_REF[veclen(n - 1, 1)];
	    REAL_REF *work_ref = new REAL_REF[max((INTEGER) 1, 2 * n - 2)];
	    REAL_REF *Dorg_ref = new REAL_REF[veclen(n, 1)];
	    REAL_REF *Eorg_ref = new REAL_REF[veclen(n, 1)];

	    COMPLEX *Z = new COMPLEX[matlen(ldz, n)];
	    REAL *D = new REAL[veclen(n, 1)];
	    REAL *E = new REAL[veclen(n - 1, 1)];
	    REAL *work = new REAL[max((INTEGER) 1, 2 * n - 2)];
	    REAL *Dorg = new REAL[veclen(n, 1)];
	    REAL *Eorg = new REAL[veclen(n - 1, 1)];
#if defined VERBOSE_TEST
	    printf("# compz %s, n:%d ldz %d\n", compz, (int) n, (int) ldz);
#endif
	    n_ref = n; ldz_ref = ldz;
	    j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(Z_ref, Z, matlen(ldz, n));
		set_random_vector(Dorg_ref, Dorg, veclen(n, 1));
		set_random_vector(Eorg_ref, Eorg, veclen(n - 1, 1));
		set_random_vector(work_ref, work, max((INTEGER) 1, 2 * n - 2));
//keep backups of D and E.
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		int iOne = 1, _n, __n;
		_n = veclen(n, 1); __n = veclen(n - 1, 1);
		dcopy_f77(&_n, Dorg_ref, &iOne, D_ref, &iOne);
		dcopy_f77(&__n, Eorg_ref, &iOne, E_ref, &iOne);
		Rcopy(veclen(n, 1), Dorg, 1, D, 1);
		Rcopy(veclen(n - 1, 1), Eorg, 1, E, 1);
#else
		Rcopy(veclen(n, 1), Dorg_ref, 1, D_ref, 1);
		Rcopy(veclen(n - 1, 1), Eorg_ref, 1, E_ref, 1);
		Rcopy(veclen(n, 1), Dorg, 1, D, 1);
		Rcopy(veclen(n - 1, 1), Eorg, 1, E, 1);
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		zsteqr_f77(compz, &n_ref, D_ref, E_ref, Z_ref, &ldz_ref, work_ref, &info_ref);
#else
		Csteqr(compz, n_ref, D_ref, E_ref, Z_ref, ldz_ref, work_ref, &info_ref);
#endif
		Csteqr(compz, n, D, E, Z, ldz, work, &info);

		if (info < 0) {
		    printf("info %d error\n", -(int) info);
		    errorflag = TRUE;
		}
		if (info_ref != info) {
		    printf("info differ! %d, %d\n", (int) info_ref, (int) info);
		    errorflag = TRUE;
		}
		if (Mlsame(compz, "I")) {
		    diff = check_eigvec <REAL_REF, COMPLEX_REF> (Dorg_ref, D_ref, Eorg_ref, Z_ref, (int) n, (int) ldz);
		}
		if (diff > EPSILON) {
		    printf("error: eigvec"); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff)
		    maxdiff = diff;
#if defined VERBOSE_TEST
		printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		if (Mlsame(compz, "I")) {
		    diff = check_eigvec <REAL, COMPLEX> (Dorg, D, Eorg, Z, (int) n, (int) ldz);
		}
		if (diff > EPSILON) {
		    printf("error: eigvec"); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff)
		    maxdiff = diff;
#if defined VERBOSE_TEST
		printf("max error: "); printnum(maxdiff); printf("\n");
#endif

		if (Mlsame(compz, "V")) {
		    /* not to be implemented */
		    /* this part can be qa'ed by Rsyev.cpp */
		}
		diff = infnorm(D_ref, D, veclen(n, 1), 1);
		if (diff > EPSILON) {
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
	    delete[]Z;
	    delete[]D;
	    delete[]E;
	    delete[]work;
	    delete[]Z_ref;
	    delete[]D_ref;
	    delete[]E_ref;
	    delete[]work_ref;
	}
	if (errorflag == TRUE) {
	    printf("*** Testing Csteqr failed ***\n");
	    exit(1);
	}
    }
}

void Csteqr_test(void)
{
    Csteqr_test2("I");
    Csteqr_test2("N");
    Csteqr_test2("V");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Csteqr start ***\n");
    Csteqr_test();
    printf("*** Testing Csteqr successful ***\n");
    return (0);
}
