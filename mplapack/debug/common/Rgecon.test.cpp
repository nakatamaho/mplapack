/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgecon.test.cpp,v 1.2 2010/08/19 01:17:55 nakatamaho Exp $
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

#define MIN_N     3
#define MAX_N     20
#define MAX_LDA   20
#define MAX_ITER  2

REAL_REF maxdiff = 0.0;

void getinvA_blas(double * A, int n, int lda)
{
    int info;
    int lwork;
/* pivot vector */
    int *ipiv = new int[veclen(n, 1)];
    lwork = -1;
    double *work = new double[1];
/* query work space */
    dgetri_f77(&n, A, &lda, ipiv, work, &lwork, &info);
    lwork = (int) (double(work[0]));
    delete[]work;
    work = new double[max(1, (int) lwork)];
/* do inversion */
    dgetrf_f77(&n, &n, A, &lda, ipiv, &info);
    dgetri_f77(&n, A, &lda, ipiv, work, &lwork, &info);
    delete[]ipiv;

    if (info == 0) return;
    if (info > 0) printf("matrix is singular\n"); exit(1);
    if (info < 0) printf("%d th argument had an illegal value\n", info); exit(1);
}

template < class REAL, class INTEGER > void getinvA(REAL * A, INTEGER n, INTEGER lda)
{
    INTEGER info;
    INTEGER lwork;
/* pivot vector */
    INTEGER *ipiv = new INTEGER[veclen(n, 1)];
    lwork = -1;
    REAL *work = new REAL[1];
/* query work space */
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
    lwork = (int) cast2double(work[0]);
    delete[]work;
    work = new REAL[max(1, (int) lwork)];
/* do inversion */
    Rgetrf(n, n, A, lda, ipiv, &info);
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
    delete[]ipiv;

    if (info == 0) return;
    if (info > 0) printf("matrix is singular\n"); exit(1);
    if (info < 0) printf("%d th argument had an illegal value\n", (int)info); exit(1);
}

void Rgecon_test2(const char *norm, int type)
{
    int errorflag = FALSE;
    int j = 0;
    REAL_REF anorm_ref;
    REAL_REF ainvnorm_ref;
    REAL_REF rcond_ref;
    REAL_REF condexact_ref;
    REAL_REF rcondexact_ref;
    REAL_REF diff;
    INTEGER_REF info_ref;
    REAL anorm;
    REAL ainvnorm;
    REAL rcond;
    REAL condexact;
    REAL rcondexact;
    INTEGER info;
    REAL_REF rtmp_ref;
    REAL rtmp;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
	    printf("#n:%d lda %d, norm %s, type %d\n", n, lda, norm, type);
#endif
	    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
	    REAL_REF *Ainv_ref = new REAL_REF[matlen(lda, n)];
	    REAL_REF *Aorg_ref = new REAL_REF[matlen(lda, n)];
	    REAL_REF *work_ref = new REAL_REF[max(1, n * 4)];
	    INTEGER_REF *iwork_ref = new INTEGER_REF[max(1, n)];
	    INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(n, 1)];

	    REAL *A = new REAL[matlen(lda, n)];
	    REAL *Ainv = new REAL[matlen(lda, n)];
	    REAL *Aorg = new REAL[matlen(lda, n)];
	    REAL *work = new REAL[max(1, n * 4)];
	    INTEGER *iwork = new INTEGER[max(1, n)];
	    INTEGER *ipiv = new INTEGER[veclen(n, 1)];

	    j = 0;
	    while (j < MAX_ITER) {
	      switch (type) {
	      case -1:
                set_random_vector(A_ref, A, matlen(lda, n));
		break;
	      case  0:
                set_hilbertmat(A_ref, A, lda, n);
		break;
	      default:
		set_random_symmmat_cond(A_ref, A, lda, n, type);
	      break;
	      }
//save A matrix
		for (int p = 0; p < n; p++) {
		    for (int q = 0; q < n; q++) {
			Ainv_ref[p + q * lda] = A_ref[p + q * lda];
			Ainv[p + q * lda] = A[p + q * lda];
		    }
		}
//get inversion of A via getri
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		getinvA_blas (Ainv_ref, n, lda);
		getinvA < REAL, INTEGER > (Ainv, n, lda);
#else
		getinvA < REAL_REF, INTEGER_REF > (Ainv_ref, n, lda);
		getinvA < REAL, INTEGER > (Ainv, n, lda);
#endif

//calculate residual
		REAL_REF *tmpmat_ref = new REAL_REF[matlen(lda, n)];
		REAL *tmpmat = new REAL[matlen(lda, n)];
	      
		for (int p = 0; p < n; p++) {
		    for (int q = 0; q < n; q++) {
		      rtmp = 0.0; rtmp_ref = 0.0;
		      for (int r = 0; r < n; r++) {
			 rtmp_ref = rtmp_ref + A_ref [ p + r * lda ] * Ainv_ref [ r + q * lda ];
			 rtmp = rtmp + A [ p + r * lda ] * Ainv [ r + q * lda ];
		      }
		    tmpmat_ref [p + q * lda] = rtmp_ref;
		    tmpmat [p + q * lda] = rtmp;
		    }
		}
	        rtmp = 0.0; rtmp_ref = 0.0;
		for (int p = 0; p < n; p++) {
		    for (int q = 0; q < n; q++) {
		      if (p!=q) {
			rtmp_ref = rtmp_ref + abs (tmpmat_ref [p + q * lda]) ;
			rtmp = rtmp + abs (tmpmat [p + q * lda]) ;
		      } else {
			rtmp_ref = rtmp_ref + abs (1.0 - tmpmat_ref [p + q * lda]) ;
			rtmp = rtmp + abs (1.0 - tmpmat [p + q * lda]) ;
		      }
		    }
		}
/*
	        printf("residual:"); printnum (rtmp); printf("\n");
	        printf("residual_ref:"); printnum (rtmp_ref); printf("\n");
*/
		delete[]tmpmat_ref;
		delete[]tmpmat;
/* calculate norm */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		anorm_ref = dlange_f77(norm, &n, &n, A_ref, &lda, work_ref);
#else
		anorm_ref = Rlange(norm, n, n, A_ref, lda, work_ref);
#endif
		anorm = Rlange(norm, n, n, A, lda, work);

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		ainvnorm_ref = dlange_f77(norm, &n, &n, Ainv_ref, &lda, work_ref);
#else
		ainvnorm_ref = Rlange(norm, n, n, Ainv_ref, lda, work_ref);
#endif
		ainvnorm = Rlange(norm, n, n, Ainv, lda, work);
/* these are very accurate condition numbers when mp is employed */
	        condexact_ref = anorm_ref * ainvnorm_ref;
	        condexact = anorm * ainvnorm;
/*
	        printf("A_ref = "); printmat (n,n,A_ref,lda); printf("\n");
	        printf("A ="); printmat (n,n,A,lda); printf("\n");
	        printf("anorm_ref:"); printnum (anorm_ref); printf("\n");
	        printf("ainvnorm_ref:"); printnum (ainvnorm_ref); printf("\n");
	        printf("anorm:"); printnum (anorm); printf("\n");
	        printf("ainvnorm:"); printnum (ainvnorm); printf("\n");
	        printf("condexact_ref:"); printnum (condexact_ref); printf("\n");
       	        printf("condexact:"); printnum (condexact); printf("\n");
*/
		diff = (condexact_ref - condexact)/condexact_ref;
		if (diff > EPSILON12) {
#if defined VERBOSE_TEST
		    printf("n:%d lda %d, norm %s\n", n, lda, norm);
		    printf("cond exact error: ");  printnum(diff); printf("\n");
	            printf("condexact_ref:"); printnum (condexact_ref); printf("\n");
       	            printf("condexact:"); printnum (condexact); printf("\n");
		    printf("\n");
#endif
		}
/* second, do LU factorization vir Rgetrf */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		dgetrf_f77(&n, &n, A_ref, &lda, ipiv_ref, &info_ref);
#else
		Rgetrf(n, n, A_ref, lda, ipiv_ref, &info_ref);
#endif
		Rgetrf(n, n, A, lda, ipiv, &info);
/* third, calculate condition number */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		dgecon_f77(norm, &n, A_ref, &lda, &anorm_ref, &rcond_ref, work_ref, iwork_ref, &info_ref);
#else
		Rgecon(norm, n, A_ref, lda, anorm_ref, &rcond_ref, work_ref, iwork_ref, &info_ref);
#endif
		Rgecon(norm, n, A, lda, anorm, &rcond, work, iwork, &info);
		if (info_ref != info) {
		    printf("info differ! %d, %d\n", (int) info_ref, (int) info);
		    errorflag = TRUE;
		}
		diff = (rcond_ref - rcond);
#if defined VERBOSE_TEST
		printf("reciprocal to cond num:"); printnum(rcond_ref); printf("\n");
		printf("reciprocal to cond num:"); printnum(rcond); printf("\n");
#endif
	        rcondexact_ref = 1.0 / condexact_ref;
	        rcondexact = 1.0 / condexact;
#if defined VERBOSE_TEST
	        printf("rcondexact_ref:"); printnum (rcondexact_ref); printf("\n");
       	        printf("rcondexact:"); printnum (rcondexact); printf("\n");
#endif
		if (diff > EPSILON) {
		    printf("n:%d lda %d, norm %s\n", n, lda, norm);
		    printf("error: ");  printnum(diff);
		    printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff) maxdiff = diff;
#if defined VERBOSE_TEST
		printf("max error: "); printnum(maxdiff);
		printf("\n");
#endif
		j++;
	    }
	    delete[]Aorg;
	    delete[]Aorg_ref;
	    delete[]Ainv;
	    delete[]Ainv_ref;
	    delete[]ipiv;
	    delete[]iwork;
	    delete[]work;
	    delete[]A;
	    delete[]ipiv_ref;
	    delete[]iwork_ref;
	    delete[]work_ref;
	    delete[]A_ref;
	}
	if (errorflag == TRUE) {
	    printf("Rgecon test failed...\n");
	    exit(1);
	}
    }
}

void Rgecon_test(void)
{
#if defined ___MPLAPACK_BUILD_WITH_GMP___ 
    Rgecon_test2("1", 100);
    Rgecon_test2("I", 100);
    Rgecon_test2("O", 100);
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___ || defined ___MPLAPACK_BUILD_WITH_QD___ || defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH__FLOAT128___
    Rgecon_test2("1", 30);
    Rgecon_test2("I", 30);
    Rgecon_test2("O", 30);
#endif

    Rgecon_test2("1", 10);
    Rgecon_test2("I", 10);
    Rgecon_test2("O", 10);

    Rgecon_test2("1", 4);
    Rgecon_test2("I", 4);
    Rgecon_test2("O", 4);

    Rgecon_test2("1", 0);
    Rgecon_test2("I", 0);
    Rgecon_test2("O", 0);

    Rgecon_test2("1", -1);
    Rgecon_test2("I", -1);
    Rgecon_test2("O", -1);
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rgecon start ***\n");
    Rgecon_test();
    printf("*** Testing Rgecon successful ***\n");
    return (0);
}
