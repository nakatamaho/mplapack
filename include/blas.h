/*
 * Copyright (c) 2008-2012
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

#ifndef _BLAS_H_
#define _BLAS_H_

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name, NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name, NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/*
 BLAS definitions are here
 Define to a macro mangling the given C identifier (in lower and upper
 case), which must not contain underscores, for linking with Fortran.
*/

#ifdef __cplusplus
extern "C" {
#endif

/* f2c.h style */
typedef int logical;

#ifdef __cplusplus
typedef std::complex<double> doublecomplex;
typedef std::complex<float> singlecomplex;
#else
typedef double _Complex doublecomplex;
typedef float _Complex singlecomplex;
#endif

#define F77_RET_I int
#define F77_RET_L logical
#define F77_RET_D double
#define F77_RET_F float
#define F77_RET_Z doublecomplex
#define F77_RET_C singlecomplex

#if defined(__APPLE__) /* Dirty... */
#define F77_FUNC(name, NAME) name ## _
#endif

#define caxpy_f77 F77_FUNC(caxpy, CAXPY)
#define ccopy_f77 F77_FUNC(ccopy, CCOPY)
#define cdotc_f77 F77_FUNC(cdotc, CDOTC)
#define cdotu_f77 F77_FUNC(cdotu, CDOTU)
#define cgbmv_f77 F77_FUNC(cgbmv, CGBMV)
#define cgemm_f77 F77_FUNC(cgemm, CGEMM)
#define cgemv_f77 F77_FUNC(cgemv, CGEMV)
#define cgerc_f77 F77_FUNC(cgerc, CGERC)
#define cgeru_f77 F77_FUNC(cgeru, CGERU)
#define chbmv_f77 F77_FUNC(chbmv, CHBMV)
#define chemm_f77 F77_FUNC(chemm, CHEMM)
#define chemv_f77 F77_FUNC(chemv, CHEMV)
#define cher_f77 F77_FUNC(cher, CHER)
#define cher2_f77 F77_FUNC(cher2, CHER2)
#define cher2k_f77 F77_FUNC(cher2k, CHER2K)
#define cherk_f77 F77_FUNC(cherk, CHERK)
#define chpmv_f77 F77_FUNC(chpmv, CHPMV)
#define chpr_f77 F77_FUNC(chpr, CHPR)
#define chpr2_f77 F77_FUNC(chpr2, CHPR2)
#define crotg_f77 F77_FUNC(crotg, CROTG)
#define cscal_f77 F77_FUNC(cscal, CSCAL)
#define csrot_f77 F77_FUNC(csrot, CSROT)
#define csscal_f77 F77_FUNC(csscal, CSSCAL)
#define cswap_f77 F77_FUNC(cswap, CSWAP)
#define csymm_f77 F77_FUNC(csymm, CSYMM)
#define csyr2k_f77 F77_FUNC(csyr2k, CSYR2K)
#define csyrk_f77 F77_FUNC(csyrk, CSYRK)
#define ctbmv_f77 F77_FUNC(ctbmv, CTBMV)
#define ctbsv_f77 F77_FUNC(ctbsv, CTBSV)
#define ctpmv_f77 F77_FUNC(ctpmv, CTPMV)
#define ctpsv_f77 F77_FUNC(ctpsv, CTPSV)
#define ctrmm_f77 F77_FUNC(ctrmm, CTRMM)
#define ctrmv_f77 F77_FUNC(ctrmv, CTRMV)
#define ctrsm_f77 F77_FUNC(ctrsm, CTRSM)
#define ctrsv_f77 F77_FUNC(ctrsv, CTRSV)
#define dasum_f77 F77_FUNC(dasum, DASUM)
#define daxpy_f77 F77_FUNC(daxpy, DAXPY)
#define dcabs1_f77 F77_FUNC(dcabs1, DCABS1)
#define dcopy_f77 F77_FUNC(dcopy, DCOPY)
#define ddot_f77 F77_FUNC(ddot, DDOT)
#define dgbmv_f77 F77_FUNC(dgbmv, DGBMV)
#define dgemm_f77 F77_FUNC(dgemm, DGEMM)
#define dgemv_f77 F77_FUNC(dgemv, DGEMV)
#define dger_f77 F77_FUNC(dger, DGER)
#define dnrm2_f77 F77_FUNC(dnrm2, DNRM2)
#define drot_f77 F77_FUNC(drot, DROT)
#define drotg_f77 F77_FUNC(drotg, DROTG)
#define drotm_f77 F77_FUNC(drotm, DROTM)
#define drotmg_f77 F77_FUNC(drotmg, DROTMG)
#define dsbmv_f77 F77_FUNC(dsbmv, DSBMV)
#define dscal_f77 F77_FUNC(dscal, DSCAL)
#define dsdot_f77 F77_FUNC(dsdot, DSDOT)
#define dspmv_f77 F77_FUNC(dspmv, DSPMV)
#define dspr_f77 F77_FUNC(dspr, DSPR)
#define dspr2_f77 F77_FUNC(dspr2, DSPR2)
#define dswap_f77 F77_FUNC(dswap, DSWAP)
#define dsymm_f77 F77_FUNC(dsymm, DSYMM)
#define dsymv_f77 F77_FUNC(dsymv, DSYMV)
#define dsyr_f77 F77_FUNC(dsyr, DSYR)
#define dsyr2_f77 F77_FUNC(dsyr2, DSYR2)
#define dsyr2k_f77 F77_FUNC(dsyr2k, DSYR2K)
#define dsyrk_f77 F77_FUNC(dsyrk, DSYRK)
#define dtbmv_f77 F77_FUNC(dtbmv, DTBMV)
#define dtbsv_f77 F77_FUNC(dtbsv, DTBSV)
#define dtpmv_f77 F77_FUNC(dtpmv, DTPMV)
#define dtpsv_f77 F77_FUNC(dtpsv, DTPSV)
#define dtrmm_f77 F77_FUNC(dtrmm, DTRMM)
#define dtrmv_f77 F77_FUNC(dtrmv, DTRMV)
#define dtrsm_f77 F77_FUNC(dtrsm, DTRSM)
#define dtrsv_f77 F77_FUNC(dtrsv, DTRSV)
#define dzasum_f77 F77_FUNC(dzasum, DZASUM)
#define dznrm2_f77 F77_FUNC(dznrm2, DZNRM2)
#define icamax_f77 F77_FUNC(icamax, ICAMAX)
#define idamax_f77 F77_FUNC(idamax, IDAMAX)
#define isamax_f77 F77_FUNC(isamax, ISAMAX)
#define izamax_f77 F77_FUNC(izamax, IZAMAX)
#define lsame_f77 F77_FUNC(lsame, LSAME)
#define sasum_f77 F77_FUNC(sasum, SASUM)
#define saxpy_f77 F77_FUNC(saxpy, SAXPY)
#define scabs1_f77 F77_FUNC(scabs1, SCABS1)
#define scasum_f77 F77_FUNC(scasum, SCASUM)
#define scnrm2_f77 F77_FUNC(scnrm2, SCNRM2)
#define scopy_f77 F77_FUNC(scopy, SCOPY)
#define sdot_f77 F77_FUNC(sdot, SDOT)
#define sdsdot_f77 F77_FUNC(sdsdot, SDSDOT)
#define sgbmv_f77 F77_FUNC(sgbmv, SGBMV)
#define sgemm_f77 F77_FUNC(sgemm, SGEMM)
#define sgemv_f77 F77_FUNC(sgemv, SGEMV)
#define sger_f77 F77_FUNC(sger, SGER)
#define snrm2_f77 F77_FUNC(snrm2, SNRM2)
#define srot_f77 F77_FUNC(srot, SROT)
#define srotg_f77 F77_FUNC(srotg, SROTG)
#define srotm_f77 F77_FUNC(srotm, SROTM)
#define srotmg_f77 F77_FUNC(srotmg, SROTMG)
#define ssbmv_f77 F77_FUNC(ssbmv, SSBMV)
#define sscal_f77 F77_FUNC(sscal, SSCAL)
#define sspmv_f77 F77_FUNC(sspmv, SSPMV)
#define sspr_f77 F77_FUNC(sspr, SSPR)
#define sspr2_f77 F77_FUNC(sspr2, SSPR2)
#define sswap_f77 F77_FUNC(sswap, SSWAP)
#define ssymm_f77 F77_FUNC(ssymm, SSYMM)
#define ssymv_f77 F77_FUNC(ssymv, SSYMV)
#define ssyr_f77 F77_FUNC(ssyr, SSYR)
#define ssyr2_f77 F77_FUNC(ssyr2, SSYR2)
#define ssyr2k_f77 F77_FUNC(ssyr2k, SSYR2K)
#define ssyrk_f77 F77_FUNC(ssyrk, SSYRK)
#define stbmv_f77 F77_FUNC(stbmv, STBMV)
#define stbsv_f77 F77_FUNC(stbsv, STBSV)
#define stpmv_f77 F77_FUNC(stpmv, STPMV)
#define stpsv_f77 F77_FUNC(stpsv, STPSV)
#define strmm_f77 F77_FUNC(strmm, STRMM)
#define strmv_f77 F77_FUNC(strmv, STRMV)
#define strsm_f77 F77_FUNC(strsm, STRSM)
#define strsv_f77 F77_FUNC(strsv, STRSV)
#define xerbla_f77 F77_FUNC(xerbla, XERBLA)
#define zaxpy_f77 F77_FUNC(zaxpy, ZAXPY)
#define zcopy_f77 F77_FUNC(zcopy, ZCOPY)
#define zdotc_f77 F77_FUNC(zdotc, ZDOTC)
#define zdotu_f77 F77_FUNC(zdotu, ZDOTU)
#define zdrot_f77 F77_FUNC(zdrot, ZDROT)
#define zdscal_f77 F77_FUNC(zdscal, ZDSCAL)
#define zgbmv_f77 F77_FUNC(zgbmv, ZGBMV)
#define zgemm_f77 F77_FUNC(zgemm, ZGEMM)
#define zgemv_f77 F77_FUNC(zgemv, ZGEMV)
#define zgerc_f77 F77_FUNC(zgerc, ZGERC)
#define zgeru_f77 F77_FUNC(zgeru, ZGERU)
#define zhbmv_f77 F77_FUNC(zhbmv, ZHBMV)
#define zhemm_f77 F77_FUNC(zhemm, ZHEMM)
#define zhemv_f77 F77_FUNC(zhemv, ZHEMV)
#define zher_f77 F77_FUNC(zher, ZHER)
#define zher2_f77 F77_FUNC(zher2, ZHER2)
#define zher2k_f77 F77_FUNC(zher2k, ZHER2K)
#define zherk_f77 F77_FUNC(zherk, ZHERK)
#define zhpmv_f77 F77_FUNC(zhpmv, ZHPMV)
#define zhpr_f77 F77_FUNC(zhpr, ZHPR)
#define zhpr2_f77 F77_FUNC(zhpr2, ZHPR2)
#define zrotg_f77 F77_FUNC(zrotg, ZROTG)
#define zscal_f77 F77_FUNC(zscal, ZSCAL)
#define zswap_f77 F77_FUNC(zswap, ZSWAP)
#define zsymm_f77 F77_FUNC(zsymm, ZSYMM)
#define zsyr2k_f77 F77_FUNC(zsyr2k, ZSYR2K)
#define zsyrk_f77 F77_FUNC(zsyrk, ZSYRK)
#define ztbmv_f77 F77_FUNC(ztbmv, ZTBMV)
#define ztbsv_f77 F77_FUNC(ztbsv, ZTBSV)
#define ztpmv_f77 F77_FUNC(ztpmv, ZTPMV)
#define ztpsv_f77 F77_FUNC(ztpsv, ZTPSV)
#define ztrmm_f77 F77_FUNC(ztrmm, ZTRMM)
#define ztrmv_f77 F77_FUNC(ztrmv, ZTRMV)
#define ztrsm_f77 F77_FUNC(ztrsm, ZTRSM)
#define ztrsv_f77 F77_FUNC(ztrsv, ZTRSV)

F77_RET_I caxpy_f77(int *n, singlecomplex *ca, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy);

F77_RET_I ccopy_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy);

F77_RET_C cdotc_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy);

F77_RET_C cdotu_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy);

F77_RET_I cgbmv_f77(const char *trans, int *m, int *n, int *kl, int *ku, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I cgemm_f77(const char *transa, const char *transb, int *m, int *n, int *k, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *beta, singlecomplex *c, int *ldc);

F77_RET_I cgemv_f77(const char *trans, int *m, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I cgerc_f77(int *m, int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *y, int *incy, singlecomplex *a, int *lda);

F77_RET_I cgeru_f77(int *m, int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *y, int *incy, singlecomplex *a, int *lda);

F77_RET_I chbmv_f77(const char *uplo, int *n, int *k, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I chemm_f77(const char *side, const char *uplo, int *m, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *beta, singlecomplex *c, int *ldc);

F77_RET_I chemv_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I cher_f77(const char *uplo, int *n, float *alpha, singlecomplex *x, int *incx, singlecomplex *a, int *lda);

F77_RET_I cher2_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *y, int *incy, singlecomplex *a, int *lda);

F77_RET_I cher2k_f77(const char *uplo, const char *trans, int *n, int *k, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *beta, singlecomplex *c, int *ldc);

F77_RET_I cherk_f77(const char *uplo, const char *trans, int *n, int *k, float *alpha, singlecomplex *a, int *lda, float *beta, singlecomplex *c, int *ldc);

F77_RET_I chpmv_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *ap, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I chpr_f77(const char *uplo, int *n, float *alpha, singlecomplex *x, int *incx, singlecomplex *ap);

F77_RET_I chpr2_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *y, int *incy, singlecomplex *ap);

F77_RET_I crotg_f77(singlecomplex *ca, singlecomplex *cb, float *c, singlecomplex *s);

F77_RET_I cscal_f77(int *n, singlecomplex *ca, singlecomplex *cx, int *incx);

F77_RET_I csrot_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy, float *c, float *s);

F77_RET_I csscal_f77(int *n, float *sa, singlecomplex *cx, int *incx);

F77_RET_I cswap_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy);

F77_RET_I csymm_f77(const char *side, const char *uplo, int *m, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *beta, singlecomplex *c, int *ldc);

F77_RET_I csyr2k_f77(const char *uplo, const char *trans, int *n, int *k, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *beta, singlecomplex *c, int *ldc);

F77_RET_I csyrk_f77(const char *uplo, const char *trans, int *n, int *k, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *beta, singlecomplex *c, int *ldc);

F77_RET_I ctbmv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *x, int *incx);

F77_RET_I ctbsv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *x, int *incx);

F77_RET_I ctpmv_f77(const char *uplo, const char *trans, const char *diag, int *n, singlecomplex *ap, singlecomplex *x, int *incx);

F77_RET_I ctpsv_f77(const char *uplo, const char *trans, const char *diag, int *n, singlecomplex *ap, singlecomplex *x, int *incx);

F77_RET_I ctrmm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb);

F77_RET_I ctrmv_f77(const char *uplo, const char *trans, const char *diag, int *n, singlecomplex *a, int *lda, singlecomplex *x, int *incx);

F77_RET_I ctrsm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *b, int *ldb);

F77_RET_I ctrsv_f77(const char *uplo, const char *trans, const char *diag, int *n, singlecomplex *a, int *lda, singlecomplex *x, int *incx);

F77_RET_D dasum_f77(int *n, double *dx, int *incx);

F77_RET_I daxpy_f77(int *n, double *da, double *dx, int *incx, double *dy, int *incy);

F77_RET_D dcabs1_f77(doublecomplex *z);

F77_RET_I dcopy_f77(int *n, double *dx, int *incx, double *dy, int *incy);

F77_RET_D ddot_f77(int *n, double *dx, int *incx, double *dy, int *incy);

F77_RET_I dgbmv_f77(const char *trans, int *m, int *n, int *kl, int *ku, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

F77_RET_I dgemm_f77(const char *transa, const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

F77_RET_I dgemv_f77(const char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

F77_RET_I dger_f77(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);

F77_RET_D dnrm2_f77(int *n, double *x, int *incx);

F77_RET_I drot_f77(int *n, double *dx, int *incx, double *dy, int *incy, double *c, double *s);

F77_RET_I drotg_f77(double *da, double *db, double *c, double *s);

F77_RET_I drotm_f77(int *n, double *dx, int *incx, double *dy, int *incy, double *dparam);

F77_RET_I drotmg_f77(double *dd1, double *dd2, double *dx1, double *dy1, double *dparam);

F77_RET_I dsbmv_f77(const char *uplo, int *n, int *k, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

F77_RET_I dscal_f77(int *n, double *da, double *dx, int *incx);

F77_RET_D dsdot_f77(int *n, float *sx, int *incx, float *sy, int *incy);

F77_RET_I dspmv_f77(const char *uplo, int *n, double *alpha, double *ap, double *x, int *incx, double *beta, double *y, int *incy);

F77_RET_I dspr_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *ap);

F77_RET_I dspr2_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *ap);

F77_RET_I dswap_f77(int *n, double *dx, int *incx, double *dy, int *incy);

F77_RET_I dsymm_f77(const char *side, const char *uplo, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

F77_RET_I dsymv_f77(const char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

F77_RET_I dsyr_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda);

F77_RET_I dsyr2_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);

F77_RET_I dsyr2k_f77(const char *uplo, const char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

F77_RET_I dsyrk_f77(const char *uplo, const char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);

F77_RET_I dtbmv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, double *a, int *lda, double *x, int *incx);

F77_RET_I dtbsv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, double *a, int *lda, double *x, int *incx);

F77_RET_I dtpmv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *ap, double *x, int *incx);

F77_RET_I dtpsv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *ap, double *x, int *incx);

F77_RET_I dtrmm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

F77_RET_I dtrmv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *a, int *lda, double *x, int *incx);

F77_RET_I dtrsm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

F77_RET_I dtrsv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *a, int *lda, double *x, int *incx);

F77_RET_D dzasum_f77(int *n, doublecomplex *zx, int *incx);

F77_RET_D dznrm2_f77(int *n, doublecomplex *x, int *incx);

F77_RET_I icamax_f77(int *n, singlecomplex *cx, int *incx);

F77_RET_I idamax_f77(int *n, double *dx, int *incx);

F77_RET_I isamax_f77(int *n, float *sx, int *incx);

F77_RET_I izamax_f77(int *n, doublecomplex *zx, int *incx);

F77_RET_L lsame_f77(const char *ca, const char *cb);

F77_RET_F sasum_f77(int *n, float *sx, int *incx);

F77_RET_I saxpy_f77(int *n, float *sa, float *sx, int *incx, float *sy, int *incy);

F77_RET_F scabs1_f77(singlecomplex *z);

F77_RET_F scasum_f77(int *n, singlecomplex *cx, int *incx);

F77_RET_F scnrm2_f77(int *n, singlecomplex *x, int *incx);

F77_RET_I scopy_f77(int *n, float *sx, int *incx, float *sy, int *incy);

F77_RET_F sdot_f77(int *n, float *sx, int *incx, float *sy, int *incy);

F77_RET_F sdsdot_f77(int *n, float *sb, float *sx, int *incx, float *sy, int *incy);

F77_RET_I sgbmv_f77(const char *trans, int *m, int *n, int *kl, int *ku, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);

F77_RET_I sgemm_f77(const char *transa, const char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

F77_RET_I sgemv_f77(const char *trans, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);

F77_RET_I sger_f77(int *m, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *a, int *lda);

F77_RET_F snrm2_f77(int *n, float *x, int *incx);

F77_RET_I srot_f77(int *n, float *sx, int *incx, float *sy, int *incy, float *c, float *s);

F77_RET_I srotg_f77(float *sa, float *sb, float *c, float *s);

F77_RET_I srotm_f77(int *n, float *sx, int *incx, float *sy, int *incy, float *sparam);

F77_RET_I srotmg_f77(float *sd1, float *sd2, float *sx1, float *sy1, float *sparam);

F77_RET_I ssbmv_f77(const char *uplo, int *n, int *k, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);

F77_RET_I sscal_f77(int *n, float *sa, float *sx, int *incx);

F77_RET_I sspmv_f77(const char *uplo, int *n, float *alpha, float *ap, float *x, int *incx, float *beta, float *y, int *incy);

F77_RET_I sspr_f77(const char *uplo, int *n, float *alpha, float *x, int *incx, float *ap);

F77_RET_I sspr2_f77(const char *uplo, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *ap);

F77_RET_I sswap_f77(int *n, float *sx, int *incx, float *sy, int *incy);

F77_RET_I ssymm_f77(const char *side, const char *uplo, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

F77_RET_I ssymv_f77(const char *uplo, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);

F77_RET_I ssyr_f77(const char *uplo, int *n, float *alpha, float *x, int *incx, float *a, int *lda);

F77_RET_I ssyr2_f77(const char *uplo, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *a, int *lda);

F77_RET_I ssyr2k_f77(const char *uplo, const char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

F77_RET_I ssyrk_f77(const char *uplo, const char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *beta, float *c, int *ldc);

F77_RET_I stbmv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, float *a, int *lda, float *x, int *incx);

F77_RET_I stbsv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, float *a, int *lda, float *x, int *incx);

F77_RET_I stpmv_f77(const char *uplo, const char *trans, const char *diag, int *n, float *ap, float *x, int *incx);

F77_RET_I stpsv_f77(const char *uplo, const char *trans, const char *diag, int *n, float *ap, float *x, int *incx);

F77_RET_I strmm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb);

F77_RET_I strmv_f77(const char *uplo, const char *trans, const char *diag, int *n, float *a, int *lda, float *x, int *incx);

F77_RET_I strsm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb);

F77_RET_I strsv_f77(const char *uplo, const char *trans, const char *diag, int *n, float *a, int *lda, float *x, int *incx);

F77_RET_I xerbla_f77(const char *srname, int *info);

F77_RET_I zaxpy_f77(int *n, doublecomplex *za, doublecomplex *zx, int *incx, doublecomplex *zy, int *incy);

F77_RET_I zcopy_f77(int *n, doublecomplex *zx, int *incx, doublecomplex *zy, int *incy);

F77_RET_Z zdotc_f77(int *n, doublecomplex *zx, int *incx, doublecomplex *zy, int *incy);

F77_RET_Z zdotu_f77(int *n, doublecomplex *zx, int *incx, doublecomplex *zy, int *incy);

F77_RET_I zdrot_f77(int *n, doublecomplex *cx, int *incx, doublecomplex *cy, int *incy, double *c, double *s);

F77_RET_I zdscal_f77(int *n, double *da, doublecomplex *zx, int *incx);

F77_RET_I zgbmv_f77(const char *trans, int *m, int *n, int *kl, int *ku, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zgemm_f77(const char *transa, const char *transb, int *m, int *n, int *k, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *beta, doublecomplex *c, int *ldc);

F77_RET_I zgemv_f77(const char *trans, int *m, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zgerc_f77(int *m, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *y, int *incy, doublecomplex *a, int *lda);

F77_RET_I zgeru_f77(int *m, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *y, int *incy, doublecomplex *a, int *lda);

F77_RET_I zhbmv_f77(const char *uplo, int *n, int *k, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zhemm_f77(const char *side, const char *uplo, int *m, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *beta, doublecomplex *c, int *ldc);

F77_RET_I zhemv_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zher_f77(const char *uplo, int *n, double *alpha, doublecomplex *x, int *incx, doublecomplex *a, int *lda);

F77_RET_I zher2_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *y, int *incy, doublecomplex *a, int *lda);

F77_RET_I zher2k_f77(const char *uplo, const char *trans, int *n, int *k, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *beta, doublecomplex *c, int *ldc);

F77_RET_I zherk_f77(const char *uplo, const char *trans, int *n, int *k, double *alpha, doublecomplex *a, int *lda, double *beta, doublecomplex *c, int *ldc);

F77_RET_I zhpmv_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *ap, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zhpr_f77(const char *uplo, int *n, double *alpha, doublecomplex *x, int *incx, doublecomplex *ap);

F77_RET_I zhpr2_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *y, int *incy, doublecomplex *ap);

F77_RET_I zrotg_f77(doublecomplex *ca, doublecomplex *cb, double *c, doublecomplex *s);

F77_RET_I zscal_f77(int *n, doublecomplex *za, doublecomplex *zx, int *incx);

F77_RET_I zswap_f77(int *n, doublecomplex *zx, int *incx, doublecomplex *zy, int *incy);

F77_RET_I zsymm_f77(const char *side, const char *uplo, int *m, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *beta, doublecomplex *c, int *ldc);

F77_RET_I zsyr2k_f77(const char *uplo, const char *trans, int *n, int *k, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *beta, doublecomplex *c, int *ldc);

F77_RET_I zsyrk_f77(const char *uplo, const char *trans, int *n, int *k, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *beta, doublecomplex *c, int *ldc);

F77_RET_I ztbmv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *x, int *incx);

F77_RET_I ztbsv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *x, int *incx);

F77_RET_I ztpmv_f77(const char *uplo, const char *trans, const char *diag, int *n, doublecomplex *ap, doublecomplex *x, int *incx);

F77_RET_I ztpsv_f77(const char *uplo, const char *trans, const char *diag, int *n, doublecomplex *ap, doublecomplex *x, int *incx);

F77_RET_I ztrmm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb);

F77_RET_I ztrmv_f77(const char *uplo, const char *trans, const char *diag, int *n, doublecomplex *a, int *lda, doublecomplex *x, int *incx);

F77_RET_I ztrsm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *b, int *ldb);

F77_RET_I ztrsv_f77(const char *uplo, const char *trans, const char *diag, int *n, doublecomplex *a, int *lda, doublecomplex *x, int *incx);

#ifdef __cplusplus
}
#endif
#endif
