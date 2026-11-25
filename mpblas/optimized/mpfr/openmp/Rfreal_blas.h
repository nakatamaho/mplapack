/*
 * Copyright (c) 2026
 *      Nakata, Maho
 *      All rights reserved.
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

// The fixed-precision MPFR BLAS on the CPU: every routine of
// MPLAPACK_FREAL_ROUTINES runs the code of the library (the reference code,
// in parallel as ../../common/blas/<name>_omp.cpp does) on copies of its
// operands in fr<PB> / fc<PB> (Rfreal_types.h) for the precisions
// PB of MPFREAL_PRECISIONS (Rgemm_freal.h), with the MPFR/MPC result bit
// for bit.  <name>_freal() (called first by <name>_omp) returns false, and
// the MPFR code then runs, for
//   - MPLAPACK_MPFR_CPU_FREAL=0, an x86-64 CPU without BMI2 and ADX,
//   - a default precision not in the list, a rounding mode other than
//     MPFR_RNDN (MPC_RNDNN), an MPC default precision other than the MPFR one,
//   - an operand of another precision, a NaN or an infinity, overlapping
//     operands, invalid arguments (the MPFR code reports them),
//   - less work than MPLAPACK_MPFR_CPU_FREAL_MIN (default 512 operations),
//   - Level 1 and 2 above MPLAPACK_MPFR_CPU_FREAL_L12_PREC bits (default 512;
//     the copies cost about as much as the operations), Level 1 with more
//     than one OpenMP thread,
//   - a result that overflows or underflows, found while running.
// Without arithmetic (Rcopy, Rswap, Ccopy, Cswap, iRamax, iCamax) or with the
// MPFR arithmetic parameters (Rnrm2, RCnrm2) the routines are not here.
// MPLAPACK_FREAL_LAPACK_ROUTINES: the same for the optimized LAPACK routines,
// run entirely in fr / fc (with a Winograd Rgemm switched on, the MPFR code).

#ifndef _MPLAPACK_RFREAL_BLAS_H_
#define _MPLAPACK_RFREAL_BLAS_H_

#include <mpblas_mpfr.h>

// X(name, (parameters), (arguments))
#define MPLAPACK_FREAL_ROUTINES(X) \
    X(Rsyrk, (const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class const beta, mpfr_class *cc, mplapackint const ldc), (uplo, trans, n, k, alpha, a, lda, beta, cc, ldc)) \
    X(Rsyr2k, (const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *cc, mplapackint const ldc), (uplo, trans, n, k, alpha, a, lda, b, ldb, beta, cc, ldc)) \
    X(Rsymm, (const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *cc, mplapackint const ldc), (side, uplo, m, n, alpha, a, lda, b, ldb, beta, cc, ldc)) \
    X(Rtrmm, (const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb), (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)) \
    X(Rtrsm, (const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb), (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)) \
    X(Rgemmtr, (const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc), (uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Rgemv, (const char *trans, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy), (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Rsymv, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy), (uplo, n, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Rtrmv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx), (uplo, trans, diag, n, a, lda, x, incx)) \
    X(Rger, (mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda), (m, n, alpha, x, incx, y, incy, a, lda)) \
    X(Rsyr, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *a, mplapackint const lda), (uplo, n, alpha, x, incx, a, lda)) \
    X(Rsyr2, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda), (uplo, n, alpha, x, incx, y, incy, a, lda)) \
    X(Rgbmv, (const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy), (trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Rsbmv, (const char *uplo, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy), (uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Rspmv, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *ap, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy), (uplo, n, alpha, ap, x, incx, beta, y, incy)) \
    X(Rspr, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *ap), (uplo, n, alpha, x, incx, ap)) \
    X(Rspr2, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *ap), (uplo, n, alpha, x, incx, y, incy, ap)) \
    X(Rtbmv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx), (uplo, trans, diag, n, k, a, lda, x, incx)) \
    X(Rtpmv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *ap, mpfr_class *x, mplapackint const incx), (uplo, trans, diag, n, ap, x, incx)) \
    X(Rtrsv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx), (uplo, trans, diag, n, a, lda, x, incx)) \
    X(Rtbsv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx), (uplo, trans, diag, n, k, a, lda, x, incx)) \
    X(Rtpsv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *ap, mpfr_class *x, mplapackint const incx), (uplo, trans, diag, n, ap, x, incx)) \
    X(Rscal, (mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx), (n, da, dx, incx)) \
    X(Rrot, (mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy, mpfr_class const c_, mpfr_class const s), (n, dx, incx, dy, incy, c_, s)) \
    X(Rasum, (mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class &result), (n, dx, incx, result)) \
    X(Cgemm, (const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc), (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Csymm, (const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc), (side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Chemm, (const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc), (side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Ctrmm, (const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb), (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)) \
    X(Ctrsm, (const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb), (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)) \
    X(Csyrk, (const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class const beta, mpc_class *c, mplapackint const ldc), (uplo, trans, n, k, alpha, a, lda, beta, c, ldc)) \
    X(Cherk, (const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpc_class *a, mplapackint const lda, mpfr_class const beta, mpc_class *c, mplapackint const ldc), (uplo, trans, n, k, alpha, a, lda, beta, c, ldc)) \
    X(Csyr2k, (const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc), (uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Cher2k, (const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpfr_class const beta, mpc_class *c, mplapackint const ldc), (uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Cgemmtr, (const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc), (uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) \
    X(Cgemv, (const char *trans, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy), (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Cgbmv, (const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy), (trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Chemv, (const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy), (uplo, n, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Chbmv, (const char *uplo, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy), (uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)) \
    X(Chpmv, (const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *ap, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy), (uplo, n, alpha, ap, x, incx, beta, y, incy)) \
    X(Cgeru, (mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda), (m, n, alpha, x, incx, y, incy, a, lda)) \
    X(Cgerc, (mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda), (m, n, alpha, x, incx, y, incy, a, lda)) \
    X(Cher, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *a, mplapackint const lda), (uplo, n, alpha, x, incx, a, lda)) \
    X(Cher2, (const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda), (uplo, n, alpha, x, incx, y, incy, a, lda)) \
    X(Chpr, (const char *uplo, mplapackint const n, mpfr_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *ap), (uplo, n, alpha, x, incx, ap)) \
    X(Chpr2, (const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *ap), (uplo, n, alpha, x, incx, y, incy, ap)) \
    X(Ctrmv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx), (uplo, trans, diag, n, a, lda, x, incx)) \
    X(Ctbmv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx), (uplo, trans, diag, n, k, a, lda, x, incx)) \
    X(Ctpmv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *ap, mpc_class *x, mplapackint const incx), (uplo, trans, diag, n, ap, x, incx)) \
    X(Ctrsv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx), (uplo, trans, diag, n, a, lda, x, incx)) \
    X(Ctbsv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx), (uplo, trans, diag, n, k, a, lda, x, incx)) \
    X(Ctpsv, (const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *ap, mpc_class *x, mplapackint const incx), (uplo, trans, diag, n, ap, x, incx)) \
    X(Caxpy, (mplapackint const n, mpc_class const za, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy), (n, za, zx, incx, zy, incy)) \
    X(Cscal, (mplapackint const n, mpc_class const za, mpc_class *zx, mplapackint const incx), (n, za, zx, incx)) \
    X(CRscal, (mplapackint const n, mpfr_class const da, mpc_class *zx, mplapackint const incx), (n, da, zx, incx)) \
    X(CRrot, (mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy, mpfr_class const c, mpfr_class const s), (n, zx, incx, zy, incy, c, s)) \
    X(Cdotc, (mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy, mpc_class &result), (n, zx, incx, zy, incy, result)) \
    X(Cdotu, (mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy, mpc_class &result), (n, zx, incx, zy, incy, result)) \
    X(RCasum, (mplapackint const n, mpc_class *zx, mplapackint const incx, mpfr_class &result), (n, zx, incx, result)) \
    X(Rdot, (mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy, mpfr_class &result), (n, dx, incx, dy, incy, result)) \
    X(Raxpy, (mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy), (n, da, dx, incx, dy, incy))


// The optimized LAPACK routines (mplapack/optimized/common), run entirely
// with fr<PB> / fc<PB> -- their BLAS calls included -- from one copy of the
// matrix: Rgesv, Rgetri, Rposv, Rpotri, ... are made of them.
#define MPLAPACK_FREAL_LAPACK_ROUTINES(X) \
    X(Rgetrf, (mplapackint const m, mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info), (m, n, a, lda, ipiv, info)) \
    X(Rgetrf2, (mplapackint const m, mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info), (m, n, a, lda, ipiv, info)) \
    X(Rpotrf, (const char *uplo, mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint &info), (uplo, n, a, lda, info)) \
    X(Rtrtri, (const char *uplo, const char *diag, mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint &info), (uplo, diag, n, a, lda, info)) \
    X(Rlauum, (const char *uplo, mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint &info), (uplo, n, a, lda, info)) \
    X(Cgetrf, (mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info), (m, n, a, lda, ipiv, info)) \
    X(Cgetrf2, (mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info), (m, n, a, lda, ipiv, info)) \
    X(Cpotrf, (const char *uplo, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint &info), (uplo, n, a, lda, info)) \
    X(Ctrtri, (const char *uplo, const char *diag, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint &info), (uplo, diag, n, a, lda, info)) \
    X(Clauum, (const char *uplo, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint &info), (uplo, n, a, lda, info)) \
    X(Rsteqr, (const char *compz, mplapackint const n, mpfr_class *d, mpfr_class *e, mpfr_class *z, mplapackint const ldz, mpfr_class *work, mplapackint &info), (compz, n, d, e, z, ldz, work, info)) \
    X(Csteqr, (const char *compz, mplapackint const n, mpfr_class *d, mpfr_class *e, mpc_class *z, mplapackint const ldz, mpfr_class *work, mplapackint &info), (compz, n, d, e, z, ldz, work, info)) \
    X(Rlasr, (const char *side, const char *pivot, const char *direct, mplapackint const m, mplapackint const n, mpfr_class *c, mpfr_class *s, mpfr_class *a, mplapackint const lda), (side, pivot, direct, m, n, c, s, a, lda)) \
    X(Clasr, (const char *side, const char *pivot, const char *direct, mplapackint const m, mplapackint const n, mpfr_class *c, mpfr_class *s, mpc_class *a, mplapackint const lda), (side, pivot, direct, m, n, c, s, a, lda))

#define MPLAPACK_FREAL_DECLARE(name, params, args) bool name##_freal params;
MPLAPACK_FREAL_ROUTINES(MPLAPACK_FREAL_DECLARE)
MPLAPACK_FREAL_LAPACK_ROUTINES(MPLAPACK_FREAL_DECLARE)
#undef MPLAPACK_FREAL_DECLARE

#endif
