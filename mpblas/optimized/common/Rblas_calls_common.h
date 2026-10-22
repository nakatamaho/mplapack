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

/*
 * From the arguments of a BLAS call to the description of the lane
 * computation of Rblas_lanes_common.h: shape, scalars and operands.  The
 * make_* functions return false, and the caller then runs the reference
 * implementation, for invalid arguments (the reference reports them), quick
 * returns, alpha == 0 (Level 2/3: the reference only scales), sizes of 0
 * and zero increments.
 *
 * An operand is X (its logical element (0,0); for a vector with a negative
 * increment, the element at the highest address), its leading dimension ld
 * (vectors: the increment), rows x cols, and the referenced region.
 */

#ifndef _MPLAPACK_RBLAS_CALLS_COMMON_H_
#define _MPLAPACK_RBLAS_CALLS_COMMON_H_

#include <algorithm>
#include <cstdlib>
#include <new>
#include "Rblas_lanes_common.h"

namespace mplapack_blas {

inline long blas_env_long(const char *name, long dflt)
{
    const char *v = std::getenv(name);
    if (v == NULL || *v == '\0')
        return dflt;
    return std::strtol(v, NULL, 10);
}

enum { REGION_ALL = 0, REGION_UPPER = 1, REGION_LOWER = 2, REGION_NODIAG = 4 };

template <class T> struct blas_operand {
    T *X;
    long ld;
    long rows, cols;
    int region;
    bool ref(long i, long j) const
    {
        if ((region & REGION_NODIAG) && i == j)
            return false;
        switch (region & 3) {
        case REGION_UPPER:
            return i <= j;
        case REGION_LOWER:
            return i >= j;
        default:
            return true;
        }
    }
};

template <class T> struct blas_call {
    blas_shape s;
    T alpha, beta;
    blas_operand<T> a, b, c; // c is the result (written); a may also be written
    bool c_read;             // the result operand is read
    bool a_written;          // Rswap, Rrot
    bool b_is_copy_of_c;     // Rtrmv: b must hold the original c
    double work;             // number of multiply-adds, for the size thresholds
};

template <class T> inline blas_operand<T> mat_operand(T *X, long ld, long rows, long cols, int region = REGION_ALL)
{
    blas_operand<T> o = {X, ld, rows, cols, region};
    return o;
}

// The vector of n elements x(0), x(inc), ..., with the BLAS convention for
// inc < 0.
template <class T> inline blas_operand<T> vec_operand(T *x, long n, long inc)
{
    blas_operand<T> o = {inc > 0 ? x : x - (n - 1) * inc, inc, 1, n, REGION_ALL};
    return o;
}

// dst = src by copy construction: an MPFR or GMP assignment would round src
// to the precision of dst, while the routines receive their scalars by value
// (copy construction keeps the precision).
template <class T> inline void blas_set_scalar(T &dst, const T &src)
{
    if (&dst == &src)
        return;
    dst.~T();
    new (&dst) T(src);
}

template <class T> inline blas_call<T> new_call(int op)
{
    blas_call<T> c;
    c.s.op = op;
    c.s.left = c.s.upper = c.s.trans = c.s.nounit = 0;
    c.s.m = c.s.n = c.s.k = 0;
    c.s.lda = c.s.ldb = c.s.ldc = 1;
    c.s.alpha_one = c.s.beta_zero = c.s.beta_one = 0;
    c.alpha = T(0.0);
    c.beta = T(0.0);
    c.a = c.b = c.c = mat_operand<T>((T *)0, 1, 0, 0);
    c.c_read = true;
    c.a_written = false;
    c.b_is_copy_of_c = false;
    c.work = 0.0;
    return c;
}

// LSAME is passed in so that this header does not depend on the type's
// mpblas header.
typedef bool (*lsame_fn)(const char *, const char *);

template <class T> bool make_syrk(blas_call<T> &c, lsame_fn ls, bool two, const char *uplo, const char *trans, long n, long k, const T &alpha, T *a, long lda, T *b, long ldb, const T &beta, T *cc, long ldc)
{
    const bool upper = ls(uplo, "U"), notr = ls(trans, "N");
    if ((!upper && !ls(uplo, "L")) || (!notr && !ls(trans, "T") && !ls(trans, "C")) || n <= 0 || k <= 0)
        return false;
    const long nrowa = notr ? n : k;
    if (lda < std::max(1L, nrowa) || (two && ldb < std::max(1L, nrowa)) || ldc < std::max(1L, n) || alpha == 0.0)
        return false;
    c = new_call<T>(two ? L3_SYR2K : L3_SYRK);
    c.s.upper = upper, c.s.trans = !notr, c.s.n = n, c.s.k = k;
    c.s.beta_zero = (beta == 0.0), c.s.beta_one = (beta == 1.0);
    blas_set_scalar(c.alpha, alpha), blas_set_scalar(c.beta, beta);
    c.a = mat_operand<T>(a, lda, nrowa, notr ? k : n);
    if (two)
        c.b = mat_operand<T>(b, ldb, nrowa, notr ? k : n);
    c.c = mat_operand<T>(cc, ldc, n, n, upper ? REGION_UPPER : REGION_LOWER);
    c.c_read = !c.s.beta_zero;
    c.work = (two ? 2.0 : 1.0) * n * n * k / 2;
    return true;
}

template <class T> bool make_symm(blas_call<T> &c, lsame_fn ls, const char *side, const char *uplo, long m, long n, const T &alpha, T *a, long lda, T *b, long ldb, const T &beta, T *cc, long ldc)
{
    const bool left = ls(side, "L"), upper = ls(uplo, "U");
    if ((!left && !ls(side, "R")) || (!upper && !ls(uplo, "L")) || m <= 0 || n <= 0)
        return false;
    const long na = left ? m : n;
    if (lda < std::max(1L, na) || ldb < std::max(1L, m) || ldc < std::max(1L, m) || alpha == 0.0)
        return false;
    c = new_call<T>(L3_SYMM);
    c.s.left = left, c.s.upper = upper, c.s.m = m, c.s.n = n;
    c.s.beta_zero = (beta == 0.0), c.s.beta_one = (beta == 1.0);
    blas_set_scalar(c.alpha, alpha), blas_set_scalar(c.beta, beta);
    c.a = mat_operand<T>(a, lda, na, na, upper ? REGION_UPPER : REGION_LOWER);
    c.b = mat_operand<T>(b, ldb, m, n);
    c.c = mat_operand<T>(cc, ldc, m, n);
    c.c_read = !c.s.beta_zero;
    c.work = (double)m * n * na;
    return true;
}

template <class T> bool make_trxm(blas_call<T> &c, lsame_fn ls, bool solve, const char *side, const char *uplo, const char *transa, const char *diag, long m, long n, const T &alpha, T *a, long lda, T *b, long ldb)
{
    const bool left = ls(side, "L"), upper = ls(uplo, "U"), nota = ls(transa, "N"), nounit = ls(diag, "N");
    if ((!left && !ls(side, "R")) || (!upper && !ls(uplo, "L")) || (!nota && !ls(transa, "T") && !ls(transa, "C")) || (!nounit && !ls(diag, "U")) || m <= 0 || n <= 0)
        return false;
    const long na = left ? m : n;
    if (lda < std::max(1L, na) || ldb < std::max(1L, m) || alpha == 0.0)
        return false;
    c = new_call<T>(solve ? L3_TRSM : L3_TRMM);
    c.s.left = left, c.s.upper = upper, c.s.trans = !nota, c.s.nounit = nounit, c.s.m = m, c.s.n = n;
    c.s.alpha_one = (alpha == 1.0);
    blas_set_scalar(c.alpha, alpha);
    c.a = mat_operand<T>(a, lda, na, na, (upper ? REGION_UPPER : REGION_LOWER) | (nounit ? 0 : REGION_NODIAG));
    c.c = mat_operand<T>(b, ldb, m, n);
    c.work = (double)m * n * na / 2;
    return true;
}

template <class T> bool make_gemv(blas_call<T> &c, lsame_fn ls, const char *trans, long m, long n, const T &alpha, T *a, long lda, T *x, long incx, const T &beta, T *y, long incy)
{
    const bool notr = ls(trans, "N");
    if ((!notr && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || lda < std::max(1L, m) || incx == 0 || incy == 0 || alpha == 0.0)
        return false;
    const long lenx = notr ? n : m, leny = notr ? m : n;
    c = new_call<T>(L2_GEMV);
    c.s.trans = !notr, c.s.m = m, c.s.n = n;
    c.s.beta_zero = (beta == 0.0), c.s.beta_one = (beta == 1.0);
    blas_set_scalar(c.alpha, alpha), blas_set_scalar(c.beta, beta);
    c.a = mat_operand<T>(a, lda, m, n);
    c.b = vec_operand<T>(x, lenx, incx);
    c.c = vec_operand<T>(y, leny, incy);
    c.c_read = !c.s.beta_zero;
    c.work = (double)m * n;
    return true;
}

template <class T> bool make_symv(blas_call<T> &c, lsame_fn ls, const char *uplo, long n, const T &alpha, T *a, long lda, T *x, long incx, const T &beta, T *y, long incy)
{
    const bool upper = ls(uplo, "U");
    if ((!upper && !ls(uplo, "L")) || n <= 0 || lda < std::max(1L, n) || incx == 0 || incy == 0 || alpha == 0.0)
        return false;
    c = new_call<T>(L2_SYMV);
    c.s.upper = upper, c.s.n = n;
    c.s.beta_zero = (beta == 0.0), c.s.beta_one = (beta == 1.0);
    blas_set_scalar(c.alpha, alpha), blas_set_scalar(c.beta, beta);
    c.a = mat_operand<T>(a, lda, n, n, upper ? REGION_UPPER : REGION_LOWER);
    c.b = vec_operand<T>(x, n, incx);
    c.c = vec_operand<T>(y, n, incy);
    c.c_read = !c.s.beta_zero;
    c.work = (double)n * n;
    return true;
}

template <class T> bool make_trmv(blas_call<T> &c, lsame_fn ls, const char *uplo, const char *trans, const char *diag, long n, T *a, long lda, T *x, long incx)
{
    const bool upper = ls(uplo, "U"), notr = ls(trans, "N"), nounit = ls(diag, "N");
    if ((!upper && !ls(uplo, "L")) || (!notr && !ls(trans, "T") && !ls(trans, "C")) || (!nounit && !ls(diag, "U")) || n <= 0 || lda < std::max(1L, n) || incx == 0)
        return false;
    c = new_call<T>(L2_TRMV);
    c.s.upper = upper, c.s.trans = !notr, c.s.nounit = nounit, c.s.n = n;
    c.a = mat_operand<T>(a, lda, n, n, (upper ? REGION_UPPER : REGION_LOWER) | (nounit ? 0 : REGION_NODIAG));
    c.c = vec_operand<T>(x, n, incx);
    c.b_is_copy_of_c = true;
    c.work = (double)n * n / 2;
    return true;
}

// Rger (sym == 0), Rsyr (sym == 1), Rsyr2 (sym == 2)
template <class T> bool make_rank(blas_call<T> &c, lsame_fn ls, int sym, const char *uplo, long m, long n, const T &alpha, T *x, long incx, T *y, long incy, T *a, long lda)
{
    bool upper = false;
    if (sym) {
        upper = ls(uplo, "U");
        if (!upper && !ls(uplo, "L"))
            return false;
        m = n;
    }
    if (m <= 0 || n <= 0 || incx == 0 || (sym != 1 && incy == 0) || lda < std::max(1L, m) || alpha == 0.0)
        return false;
    c = new_call<T>(sym == 0 ? L2_GER : sym == 1 ? L2_SYR : L2_SYR2);
    c.s.upper = upper, c.s.m = m, c.s.n = n;
    blas_set_scalar(c.alpha, alpha);
    c.a = vec_operand<T>(x, m, incx);
    if (sym != 1)
        c.b = vec_operand<T>(y, n, incy);
    c.c = mat_operand<T>(a, lda, m, n, sym ? (upper ? REGION_UPPER : REGION_LOWER) : REGION_ALL);
    c.work = (double)m * n;
    return true;
}

// Rscal (alpha = da), Raxpy (alpha = da), Rcopy, Rswap, Rrot (alpha = c, beta = s)
template <class T> bool make_level1(blas_call<T> &c, int op, long n, const T &alpha, const T &beta, T *x, long incx, T *y, long incy)
{
    if (n <= 0 || incx == 0 || (op != L1_SCAL && incy == 0))
        return false;
    if (op == L1_SCAL && (incx <= 0 || alpha == 1.0)) // Rscal returns for these
        return false;
    if (op == L1_AXPY && alpha == 0.0)
        return false;
    c = new_call<T>(op);
    c.s.n = n;
    blas_set_scalar(c.alpha, alpha), blas_set_scalar(c.beta, beta);
    if (op == L1_SCAL) {
        c.c = vec_operand<T>(x, n, incx);
    } else {
        c.a = vec_operand<T>(x, n, incx);
        c.c = vec_operand<T>(y, n, incy);
        c.c_read = op != L1_COPY;
        c.a_written = op == L1_SWAP || op == L1_ROT;
        // the parallel lanes need x and y not to overlap
        const T *x0 = std::min(c.a.X, c.a.X + (n - 1) * incx), *x1 = std::max(c.a.X, c.a.X + (n - 1) * incx);
        const T *y0 = std::min(c.c.X, c.c.X + (n - 1) * incy), *y1 = std::max(c.c.X, c.c.X + (n - 1) * incy);
        if (!(x1 < y0 || y1 < x0))
            return false;
    }
    c.work = (double)n;
    return true;
}

} // namespace mplapack_blas

#endif
