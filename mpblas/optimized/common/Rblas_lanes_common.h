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
 * Level 1, 2 and 3 BLAS written once against an arithmetic engine, for the
 * CPU (OpenMP, any MPLAPACK type) and for CUDA (dd, td, qd, MPFR).
 *
 * Every routine performs, for every element of its result, the operations
 * of the reference BLAS (mpblas/reference/R*.cpp) in the same order, so that
 * with an engine whose operations are those of the C++ type (or reproduce
 * them bit for bit) the result is the reference result exactly.  The work
 * is split into independent lanes (see blas_lanes()):
 *
 *   one element      Rsyrk, Rsyr2k (triangle of C), Rger, Rsyr, Rsyr2
 *                    (A), Rgemv, Rsymv, Rtrmv (y / x), Rscal, Raxpy,
 *                    Rcopy, Rswap, Rrot
 *   one column       Rsymm, Rtrmm, Rtrsm with side = L
 *   one row          Rsymm, Rtrmm, Rtrsm with side = R
 *
 * Statements are decomposed into single operations in the order the C++
 * types evaluate them: "c += a*t1 + b*t2" is c = c + ((a*t1) + (b*t2)), and
 * "x = p*q + r*s + u*v" is x = ((p*q) + (r*s)) + (u*v); the operands of
 * every multiplication keep their order (qd_real multiplication is not
 * commutative).
 *
 * Engine interface (the functions are MPLAPACK_BLAS_HD for CUDA engines):
 *   typedef ... V;  struct M { ... };          value / matrix handle
 *   V reg(int slot);                            a register (slot < 8)
 *   at(const M &, long i, long j)              read-only value (V or const V &)
 *   void put(M &, long i, long j, const V &);
 *   void set(V &, const V &); void zero(V &); void one(V &);
 *   void add(V &, const V &, const V &);  sub, mul, div
 *   bool iszero(const V &); bool isone(const V &);
 *   static const bool has_div;
 * A vector is a 1 x n matrix whose leading dimension is the increment
 * (element i at column i); M of a vector points at its logical element 0.
 */

#ifndef _MPLAPACK_RBLAS_LANES_COMMON_H_
#define _MPLAPACK_RBLAS_LANES_COMMON_H_

#if defined(__CUDACC__) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#define MPLAPACK_BLAS_HD __host__ __device__
#else
#define MPLAPACK_BLAS_HD
#endif

namespace mplapack_blas {

enum {
    L3_SYRK = 1, L3_SYR2K = 2, L3_SYMM = 3, L3_TRMM = 4, L3_TRSM = 5,
    L2_GEMV = 11, L2_SYMV = 12, L2_TRMV = 13, L2_GER = 14, L2_SYR = 15, L2_SYR2 = 16,
    L1_SCAL = 21, L1_AXPY = 22, L1_COPY = 23, L1_SWAP = 24, L1_ROT = 25
};

// Operands of each routine:
//   Level 3   a = A, b = B, c = C (Rtrmm/Rtrsm: B); al = alpha, be = beta
//   GEMV      a = A, b = x, c = y        SYMV  a = A, b = x, c = y
//   TRMV      a = A, b = copy of x (read), c = x
//   GER       a = x, b = y, c = A        SYR   a = x, c = A      SYR2  a = x, b = y, c = A
//   SCAL      c = x (al = alpha)         AXPY  a = x, c = y (al = alpha)
//   COPY      a = x, c = y               SWAP  a = x, c = y (both written)
//   ROT       a = x, c = y (both written), al = c, be = s
struct blas_shape {
    int op;
    int left;     // side == "L"
    int upper;    // uplo == "U"
    int trans;    // trans / transa != "N"
    int nounit;   // diag == "N"
    long m, n, k; // syrk/syr2k: n, k; Level 2: m, n (square: n); Level 1: n
    long lda, ldb, ldc;
    int alpha_one, beta_zero, beta_one;
};

// Number of lanes of a call.
MPLAPACK_BLAS_HD inline long blas_lanes(const blas_shape &s)
{
    switch (s.op) {
    case L3_SYRK:
    case L3_SYR2K:
    case L2_SYR:
    case L2_SYR2:
        return s.n * s.n;
    case L3_SYMM:
    case L3_TRMM:
    case L3_TRSM:
        return s.left ? s.n : s.m;
    case L2_GEMV:
        return s.trans ? s.n : s.m;
    case L2_GER:
        return s.m * s.n;
    default: // SYMV, TRMV, Level 1
        return s.n;
    }
}

template <class E> struct blas_mats {
    typename E::M al, be, a, b, c; // alpha and beta are 1 x 1
};

// ---- Rsyrk / Rsyr2k: element (i, j) of the referenced triangle ----

template <class E> MPLAPACK_BLAS_HD void l3_syrk(E &e, const blas_shape &s, blas_mats<E> &x, long i, long j)
{
    typedef typename E::V V;
    if (s.upper ? i > j : i < j)
        return;
    const V al = e.at(x.al, 0, 0), be = e.at(x.be, 0, 0);
    V c = e.reg(0), p = e.reg(1), t = e.reg(2);
    if (!s.trans) {
        if (s.beta_zero)
            e.zero(c);
        else if (!s.beta_one)
            e.mul(c, be, e.at(x.c, i, j));
        else
            e.set(c, e.at(x.c, i, j));
        for (long l = 0; l < s.k; l++) {
            if (!e.iszero(e.at(x.a, j, l))) {
                e.mul(t, al, e.at(x.a, j, l));   // temp = alpha * a(j,l)
                e.mul(p, t, e.at(x.a, i, l));    // c(i,j) += temp * a(i,l)
                e.add(c, c, p);
            }
        }
    } else {
        e.zero(t);
        for (long l = 0; l < s.k; l++) {
            e.mul(p, e.at(x.a, l, i), e.at(x.a, l, j)); // temp += a(l,i) * a(l,j)
            e.add(t, t, p);
        }
        e.mul(c, al, t);                         // c = alpha*temp [+ beta*c]
        if (!s.beta_zero) {
            e.mul(p, be, e.at(x.c, i, j));
            e.add(c, c, p);
        }
    }
    e.put(x.c, i, j, c);
}

template <class E> MPLAPACK_BLAS_HD void l3_syr2k(E &e, const blas_shape &s, blas_mats<E> &x, long i, long j)
{
    typedef typename E::V V;
    if (s.upper ? i > j : i < j)
        return;
    const V al = e.at(x.al, 0, 0), be = e.at(x.be, 0, 0);
    V c = e.reg(0), p = e.reg(1), q = e.reg(2), t1 = e.reg(3), t2 = e.reg(4);
    if (!s.trans) {
        if (s.beta_zero)
            e.zero(c);
        else if (!s.beta_one)
            e.mul(c, be, e.at(x.c, i, j));
        else
            e.set(c, e.at(x.c, i, j));
        for (long l = 0; l < s.k; l++) {
            if (!e.iszero(e.at(x.a, j, l)) || !e.iszero(e.at(x.b, j, l))) {
                e.mul(t1, al, e.at(x.b, j, l)); // temp1 = alpha * b(j,l)
                e.mul(t2, al, e.at(x.a, j, l)); // temp2 = alpha * a(j,l)
                e.mul(p, e.at(x.a, i, l), t1);  // c += a(i,l)*temp1 + b(i,l)*temp2
                e.mul(q, e.at(x.b, i, l), t2);
                e.add(p, p, q);
                e.add(c, c, p);
            }
        }
    } else {
        e.zero(t1);
        e.zero(t2);
        for (long l = 0; l < s.k; l++) {
            e.mul(p, e.at(x.a, l, i), e.at(x.b, l, j)); // temp1 += a(l,i) * b(l,j)
            e.add(t1, t1, p);
            e.mul(p, e.at(x.b, l, i), e.at(x.a, l, j)); // temp2 += b(l,i) * a(l,j)
            e.add(t2, t2, p);
        }
        if (s.beta_zero) {
            e.mul(c, al, t1); // alpha*temp1 + alpha*temp2
            e.mul(p, al, t2);
            e.add(c, c, p);
        } else {
            e.mul(c, be, e.at(x.c, i, j)); // beta*c + alpha*temp1 + alpha*temp2
            e.mul(p, al, t1);
            e.add(c, c, p);
            e.mul(p, al, t2);
            e.add(c, c, p);
        }
    }
    e.put(x.c, i, j, c);
}

// ---- Rsymm: column j (side L) or row i (side R) of C ----

template <class E> MPLAPACK_BLAS_HD void l3_symm(E &e, const blas_shape &s, blas_mats<E> &x, long lane)
{
    typedef typename E::V V;
    const V al = e.at(x.al, 0, 0), be = e.at(x.be, 0, 0);
    V c = e.reg(0), p = e.reg(1), t1 = e.reg(2), t2 = e.reg(3), y = e.reg(4);
    if (s.left) {
        const long j = lane;
        for (long ii = 0; ii < s.m; ii++) {
            const long i = s.upper ? ii : s.m - 1 - ii;
            e.mul(t1, al, e.at(x.b, i, j)); // temp1 = alpha * b(i,j)
            e.zero(t2);
            const long k0 = s.upper ? 0 : i + 1, k1 = s.upper ? i : s.m;
            for (long k = k0; k < k1; k++) {
                e.set(y, e.at(x.c, k, j)); // c(k,j) += temp1 * a(k,i)
                e.mul(p, t1, e.at(x.a, k, i));
                e.add(y, y, p);
                e.put(x.c, k, j, y);
                e.mul(p, e.at(x.b, k, j), e.at(x.a, k, i)); // temp2 += b(k,j) * a(k,i)
                e.add(t2, t2, p);
            }
            if (s.beta_zero) {
                e.mul(c, t1, e.at(x.a, i, i)); // temp1*a(i,i) + alpha*temp2
            } else {
                e.mul(c, be, e.at(x.c, i, j)); // beta*c + temp1*a(i,i) + alpha*temp2
                e.mul(p, t1, e.at(x.a, i, i));
                e.add(c, c, p);
            }
            e.mul(p, al, t2);
            e.add(c, c, p);
            e.put(x.c, i, j, c);
        }
    } else {
        const long i = lane;
        for (long j = 0; j < s.n; j++) {
            e.mul(t1, al, e.at(x.a, j, j)); // temp1 = alpha * a(j,j)
            if (s.beta_zero) {
                e.mul(c, t1, e.at(x.b, i, j));
            } else {
                e.mul(c, be, e.at(x.c, i, j)); // beta*c + temp1*b(i,j)
                e.mul(p, t1, e.at(x.b, i, j));
                e.add(c, c, p);
            }
            for (long k = 0; k < j; k++) {
                e.mul(t1, al, s.upper ? e.at(x.a, k, j) : e.at(x.a, j, k));
                e.mul(p, t1, e.at(x.b, i, k)); // c(i,j) += temp1 * b(i,k)
                e.add(c, c, p);
            }
            for (long k = j + 1; k < s.n; k++) {
                e.mul(t1, al, s.upper ? e.at(x.a, j, k) : e.at(x.a, k, j));
                e.mul(p, t1, e.at(x.b, i, k));
                e.add(c, c, p);
            }
            e.put(x.c, i, j, c);
        }
    }
}

// b(i,j) += t * w
template <class E> MPLAPACK_BLAS_HD void l3_axpy1(E &e, typename E::M &B, long i, long j, const typename E::V &t, const typename E::V &w, typename E::V &y, typename E::V &p)
{
    e.set(y, e.at(B, i, j));
    e.mul(p, t, w);
    e.add(y, y, p);
    e.put(B, i, j, y);
}

// b(i,j) = b(i,j) - t * w
template <class E> MPLAPACK_BLAS_HD void l3_axmy1(E &e, typename E::M &B, long i, long j, const typename E::V &t, const typename E::V &w, typename E::V &y, typename E::V &p)
{
    e.set(y, e.at(B, i, j));
    e.mul(p, t, w);
    e.sub(y, y, p);
    e.put(B, i, j, y);
}

// b(i,j) = t * b(i,j)
template <class E> MPLAPACK_BLAS_HD void l3_scale1(E &e, typename E::M &B, long i, long j, const typename E::V &t, typename E::V &y)
{
    e.set(y, e.at(B, i, j));
    e.mul(y, t, y);
    e.put(B, i, j, y);
}

// ---- Rtrmm: B := alpha*op(A)*B (side L, column lane) or alpha*B*op(A) (side R, row lane) ----
// B is x.c.

template <class E> MPLAPACK_BLAS_HD void l3_trmm(E &e, const blas_shape &s, blas_mats<E> &x, long lane)
{
    typedef typename E::V V;
    const V al = e.at(x.al, 0, 0);
    V t = e.reg(0), p = e.reg(1), y = e.reg(2);
    const long m = s.m, n = s.n;
    if (s.left) {
        const long j = lane;
        if (!s.trans) {
            for (long kk = 0; kk < m; kk++) {
                const long k = s.upper ? kk : m - 1 - kk;
                if (e.iszero(e.at(x.c, k, j)))
                    continue;
                e.mul(t, al, e.at(x.c, k, j)); // temp = alpha * b(k,j)
                if (s.upper) {
                    for (long i = 0; i < k; i++)
                        l3_axpy1(e, x.c, i, j, t, e.at(x.a, i, k), y, p);
                    if (s.nounit)
                        e.mul(t, t, e.at(x.a, k, k));
                    e.put(x.c, k, j, t);
                } else {
                    if (s.nounit) {
                        e.mul(y, t, e.at(x.a, k, k)); // b(k,j) = temp; b(k,j) = b(k,j)*a(k,k)
                        e.put(x.c, k, j, y);
                    } else {
                        e.put(x.c, k, j, t);
                    }
                    for (long i = k + 1; i < m; i++)
                        l3_axpy1(e, x.c, i, j, t, e.at(x.a, i, k), y, p);
                }
            }
        } else {
            for (long ii = 0; ii < m; ii++) {
                const long i = s.upper ? m - 1 - ii : ii;
                e.set(t, e.at(x.c, i, j)); // temp = b(i,j)
                if (s.nounit)
                    e.mul(t, t, e.at(x.a, i, i));
                const long k0 = s.upper ? 0 : i + 1, k1 = s.upper ? i : m;
                for (long k = k0; k < k1; k++) {
                    e.mul(p, e.at(x.a, k, i), e.at(x.c, k, j)); // temp += a(k,i) * b(k,j)
                    e.add(t, t, p);
                }
                e.mul(y, al, t); // b(i,j) = alpha * temp
                e.put(x.c, i, j, y);
            }
        }
    } else {
        const long i = lane;
        if (!s.trans) {
            for (long jj = 0; jj < n; jj++) {
                const long j = s.upper ? n - 1 - jj : jj;
                e.set(t, al); // temp = alpha [* a(j,j)]
                if (s.nounit)
                    e.mul(t, t, e.at(x.a, j, j));
                l3_scale1(e, x.c, i, j, t, y);
                const long k0 = s.upper ? 0 : j + 1, k1 = s.upper ? j : n;
                for (long k = k0; k < k1; k++) {
                    if (!e.iszero(e.at(x.a, k, j))) {
                        e.mul(t, al, e.at(x.a, k, j)); // temp = alpha * a(k,j)
                        l3_axpy1(e, x.c, i, j, t, e.at(x.c, i, k), y, p);
                    }
                }
            }
        } else {
            for (long kk = 0; kk < n; kk++) {
                const long k = s.upper ? kk : n - 1 - kk;
                const long j0 = s.upper ? 0 : k + 1, j1 = s.upper ? k : n;
                for (long j = j0; j < j1; j++) {
                    if (!e.iszero(e.at(x.a, j, k))) {
                        e.mul(t, al, e.at(x.a, j, k)); // temp = alpha * a(j,k)
                        l3_axpy1(e, x.c, i, j, t, e.at(x.c, i, k), y, p);
                    }
                }
                e.set(t, al);
                if (s.nounit)
                    e.mul(t, t, e.at(x.a, k, k));
                if (!e.isone(t))
                    l3_scale1(e, x.c, i, k, t, y);
            }
        }
    }
}

// ---- Rtrsm: solves op(A)*X = alpha*B (side L) or X*op(A) = alpha*B (side R) ----

template <class E> MPLAPACK_BLAS_HD void l3_trsm(E &e, const blas_shape &s, blas_mats<E> &x, long lane)
{
    typedef typename E::V V;
    const V al = e.at(x.al, 0, 0);
    V t = e.reg(0), p = e.reg(1), y = e.reg(2), bk = e.reg(3), o = e.reg(4);
    const long m = s.m, n = s.n;
    e.one(o);
    if (s.left) {
        const long j = lane;
        if (!s.trans) {
            if (!s.alpha_one)
                for (long i = 0; i < m; i++)
                    l3_scale1(e, x.c, i, j, al, y);
            for (long kk = 0; kk < m; kk++) {
                const long k = s.upper ? m - 1 - kk : kk;
                if (e.iszero(e.at(x.c, k, j)))
                    continue;
                if (s.nounit) {
                    e.div(y, e.at(x.c, k, j), e.at(x.a, k, k)); // b(k,j) = b(k,j) / a(k,k)
                    e.put(x.c, k, j, y);
                }
                e.set(bk, e.at(x.c, k, j));
                const long i0 = s.upper ? 0 : k + 1, i1 = s.upper ? k : m;
                for (long i = i0; i < i1; i++)
                    l3_axmy1(e, x.c, i, j, bk, e.at(x.a, i, k), y, p); // b(i,j) -= b(k,j)*a(i,k)
            }
        } else {
            for (long ii = 0; ii < m; ii++) {
                const long i = s.upper ? ii : m - 1 - ii;
                e.mul(t, al, e.at(x.c, i, j)); // temp = alpha * b(i,j)
                const long k0 = s.upper ? 0 : i + 1, k1 = s.upper ? i : m;
                for (long k = k0; k < k1; k++) {
                    e.mul(p, e.at(x.a, k, i), e.at(x.c, k, j)); // temp -= a(k,i) * b(k,j)
                    e.sub(t, t, p);
                }
                if (s.nounit)
                    e.div(t, t, e.at(x.a, i, i));
                e.put(x.c, i, j, t);
            }
        }
    } else {
        const long i = lane;
        if (!s.trans) {
            for (long jj = 0; jj < n; jj++) {
                const long j = s.upper ? jj : n - 1 - jj;
                if (!s.alpha_one)
                    l3_scale1(e, x.c, i, j, al, y);
                const long k0 = s.upper ? 0 : j + 1, k1 = s.upper ? j : n;
                for (long k = k0; k < k1; k++)
                    if (!e.iszero(e.at(x.a, k, j)))
                        l3_axmy1(e, x.c, i, j, e.at(x.a, k, j), e.at(x.c, i, k), y, p); // b(i,j) -= a(k,j)*b(i,k)
                if (s.nounit) {
                    e.div(t, o, e.at(x.a, j, j)); // temp = one / a(j,j)
                    l3_scale1(e, x.c, i, j, t, y);
                }
            }
        } else {
            for (long kk = 0; kk < n; kk++) {
                const long k = s.upper ? n - 1 - kk : kk;
                if (s.nounit) {
                    e.div(t, o, e.at(x.a, k, k)); // temp = one / a(k,k)
                    l3_scale1(e, x.c, i, k, t, y);
                }
                const long j0 = s.upper ? 0 : k + 1, j1 = s.upper ? k : n;
                for (long j = j0; j < j1; j++)
                    if (!e.iszero(e.at(x.a, j, k)))
                        l3_axmy1(e, x.c, i, j, e.at(x.a, j, k), e.at(x.c, i, k), y, p); // b(i,j) -= temp*b(i,k)
                if (!s.alpha_one)
                    l3_scale1(e, x.c, i, k, al, y);
            }
        }
    }
}

// ---- Level 2 ----

// y(i) := beta*y(i) as the first loop of Rgemv / Rsymv (y unchanged when beta == 1)
template <class E> MPLAPACK_BLAS_HD void blas_beta_y(E &e, const blas_shape &s, blas_mats<E> &x, long i, typename E::V &y)
{
    if (s.beta_one)
        e.set(y, e.at(x.c, 0, i));
    else if (s.beta_zero)
        e.zero(y);
    else
        e.mul(y, e.at(x.be, 0, 0), e.at(x.c, 0, i));
}

template <class E> MPLAPACK_BLAS_HD void l2_gemv(E &e, const blas_shape &s, blas_mats<E> &x, long lane)
{
    typedef typename E::V V;
    const V al = e.at(x.al, 0, 0);
    V y = e.reg(0), p = e.reg(1), t = e.reg(2);
    blas_beta_y(e, s, x, lane, y);
    if (!s.trans) {
        const long i = lane;
        for (long j = 0; j < s.n; j++) {
            e.mul(t, al, e.at(x.b, 0, j));    // temp = alpha * x(j)
            e.mul(p, t, e.at(x.a, i, j));     // y(i) += temp * a(i,j)
            e.add(y, y, p);
        }
    } else {
        const long j = lane;
        e.zero(t);
        for (long i = 0; i < s.m; i++) {
            e.mul(p, e.at(x.a, i, j), e.at(x.b, 0, i)); // temp += a(i,j) * x(i)
            e.add(t, t, p);
        }
        e.mul(p, al, t); // y(j) += alpha * temp
        e.add(y, y, p);
    }
    e.put(x.c, 0, lane, y);
}

template <class E> MPLAPACK_BLAS_HD void l2_symv(E &e, const blas_shape &s, blas_mats<E> &x, long i)
{
    typedef typename E::V V;
    const V al = e.at(x.al, 0, 0);
    V y = e.reg(0), p = e.reg(1), q = e.reg(2), t1 = e.reg(3), t2 = e.reg(4);
    blas_beta_y(e, s, x, i, y);
    if (s.upper) {
        // j = i: y(i) += temp1*a(i,i) + alpha*temp2, temp2 = sum_{k<i} a(k,i)*x(k)
        e.mul(t1, al, e.at(x.b, 0, i));
        e.zero(t2);
        for (long k = 0; k < i; k++) {
            e.mul(p, e.at(x.a, k, i), e.at(x.b, 0, k));
            e.add(t2, t2, p);
        }
        e.mul(p, t1, e.at(x.a, i, i));
        e.mul(q, al, t2);
        e.add(p, p, q);
        e.add(y, y, p);
        // j > i: y(i) += temp1(j) * a(i,j)
        for (long j = i + 1; j < s.n; j++) {
            e.mul(t1, al, e.at(x.b, 0, j));
            e.mul(p, t1, e.at(x.a, i, j));
            e.add(y, y, p);
        }
    } else {
        // j < i: y(i) += temp1(j) * a(i,j)
        for (long j = 0; j < i; j++) {
            e.mul(t1, al, e.at(x.b, 0, j));
            e.mul(p, t1, e.at(x.a, i, j));
            e.add(y, y, p);
        }
        // j = i: y(i) += temp1*a(i,i); ...; y(i) += alpha*temp2, temp2 = sum_{k>i} a(k,i)*x(k)
        e.mul(t1, al, e.at(x.b, 0, i));
        e.mul(p, t1, e.at(x.a, i, i));
        e.add(y, y, p);
        e.zero(t2);
        for (long k = i + 1; k < s.n; k++) {
            e.mul(p, e.at(x.a, k, i), e.at(x.b, 0, k));
            e.add(t2, t2, p);
        }
        e.mul(p, al, t2);
        e.add(y, y, p);
    }
    e.put(x.c, 0, i, y);
}

// x := op(A) x; x.b is the original x, x.c receives the result.
template <class E> MPLAPACK_BLAS_HD void l2_trmv(E &e, const blas_shape &s, blas_mats<E> &x, long i)
{
    typedef typename E::V V;
    V y = e.reg(0), p = e.reg(1);
    const long n = s.n;
    if (!s.trans) {
        // iteration j = i (only when x(i) != 0): x(i) = x(i) * a(i,i);
        // then x(i) += x(j) * a(i,j) for the later iterations j that update i
        e.set(y, e.at(x.b, 0, i));
        if (s.nounit && !e.iszero(e.at(x.b, 0, i)))
            e.mul(y, y, e.at(x.a, i, i));
        if (s.upper) {
            for (long j = i + 1; j < n; j++)
                if (!e.iszero(e.at(x.b, 0, j))) {
                    e.mul(p, e.at(x.b, 0, j), e.at(x.a, i, j)); // x(i) += temp * a(i,j)
                    e.add(y, y, p);
                }
        } else {
            for (long j = i - 1; j >= 0; j--)
                if (!e.iszero(e.at(x.b, 0, j))) {
                    e.mul(p, e.at(x.b, 0, j), e.at(x.a, i, j));
                    e.add(y, y, p);
                }
        }
    } else {
        const long j = i;
        e.set(y, e.at(x.b, 0, j)); // temp = x(j) [* a(j,j)]
        if (s.nounit)
            e.mul(y, y, e.at(x.a, j, j));
        if (s.upper) {
            for (long k = j - 1; k >= 0; k--) {
                e.mul(p, e.at(x.a, k, j), e.at(x.b, 0, k)); // temp += a(k,j) * x(k)
                e.add(y, y, p);
            }
        } else {
            for (long k = j + 1; k < n; k++) {
                e.mul(p, e.at(x.a, k, j), e.at(x.b, 0, k));
                e.add(y, y, p);
            }
        }
    }
    e.put(x.c, 0, i, y);
}

// Rger (all elements), Rsyr / Rsyr2 (triangle): element (i, j) of A
template <class E> MPLAPACK_BLAS_HD void l2_rank(E &e, const blas_shape &s, blas_mats<E> &x, long i, long j)
{
    typedef typename E::V V;
    if (s.op != L2_GER && (s.upper ? i > j : i < j))
        return;
    const V al = e.at(x.al, 0, 0);
    V y = e.reg(0), p = e.reg(1), q = e.reg(2), t1 = e.reg(3), t2 = e.reg(4);
    if (s.op == L2_SYR2) {
        if (e.iszero(e.at(x.a, 0, j)) && e.iszero(e.at(x.b, 0, j)))
            return;
        e.mul(t1, al, e.at(x.b, 0, j)); // temp1 = alpha * y(j)
        e.mul(t2, al, e.at(x.a, 0, j)); // temp2 = alpha * x(j)
        e.mul(p, e.at(x.a, 0, i), t1);  // a(i,j) += x(i)*temp1 + y(i)*temp2
        e.mul(q, e.at(x.b, 0, i), t2);
        e.add(p, p, q);
    } else {
        const typename E::M &w = s.op == L2_GER ? x.b : x.a; // Rger: y(j); Rsyr: x(j)
        if (e.iszero(e.at(w, 0, j)))
            return;
        e.mul(t1, al, e.at(w, 0, j));  // temp = alpha * y(j) (Rsyr: x(j))
        e.mul(p, e.at(x.a, 0, i), t1); // a(i,j) += x(i) * temp
    }
    e.set(y, e.at(x.c, i, j));
    e.add(y, y, p);
    e.put(x.c, i, j, y);
}

// ---- Level 1: element i ----

template <class E> MPLAPACK_BLAS_HD void l1_element(E &e, const blas_shape &s, blas_mats<E> &x, long i)
{
    typedef typename E::V V;
    V y = e.reg(0), p = e.reg(1), q = e.reg(2);
    switch (s.op) {
    case L1_SCAL: // dx(i) = da * dx(i)
        e.mul(y, e.at(x.al, 0, 0), e.at(x.c, 0, i));
        e.put(x.c, 0, i, y);
        break;
    case L1_AXPY: // dy(i) += da * dx(i)
        e.mul(p, e.at(x.al, 0, 0), e.at(x.a, 0, i));
        e.set(y, e.at(x.c, 0, i));
        e.add(y, y, p);
        e.put(x.c, 0, i, y);
        break;
    case L1_COPY:
        e.set(y, e.at(x.a, 0, i));
        e.put(x.c, 0, i, y);
        break;
    case L1_SWAP:
        e.set(y, e.at(x.a, 0, i));
        e.set(p, e.at(x.c, 0, i));
        e.put(x.a, 0, i, p);
        e.put(x.c, 0, i, y);
        break;
    case L1_ROT: // dtemp = c*dx + s*dy; dy = c*dy - s*dx; dx = dtemp
        e.mul(y, e.at(x.al, 0, 0), e.at(x.a, 0, i));
        e.mul(p, e.at(x.be, 0, 0), e.at(x.c, 0, i));
        e.add(y, y, p);
        e.mul(p, e.at(x.al, 0, 0), e.at(x.c, 0, i));
        e.mul(q, e.at(x.be, 0, 0), e.at(x.a, 0, i));
        e.sub(p, p, q);
        e.put(x.c, 0, i, p);
        e.put(x.a, 0, i, y);
        break;
    default:
        break;
    }
}

template <class E> MPLAPACK_BLAS_HD void blas_lane(E &e, const blas_shape &s, blas_mats<E> &x, long lane)
{
    switch (s.op) {
    case L3_SYRK:
        l3_syrk(e, s, x, lane % s.n, lane / s.n);
        break;
    case L3_SYR2K:
        l3_syr2k(e, s, x, lane % s.n, lane / s.n);
        break;
    case L3_SYMM:
        l3_symm(e, s, x, lane);
        break;
    case L3_TRMM:
        l3_trmm(e, s, x, lane);
        break;
    case L3_TRSM:
        if (E::has_div)
            l3_trsm(e, s, x, lane);
        break;
    case L2_GEMV:
        l2_gemv(e, s, x, lane);
        break;
    case L2_SYMV:
        l2_symv(e, s, x, lane);
        break;
    case L2_TRMV:
        l2_trmv(e, s, x, lane);
        break;
    case L2_GER:
        l2_rank(e, s, x, lane % s.m, lane / s.m);
        break;
    case L2_SYR:
    case L2_SYR2:
        l2_rank(e, s, x, lane % s.n, lane / s.n);
        break;
    default:
        l1_element(e, s, x, lane);
        break;
    }
}

} // namespace mplapack_blas

#endif
