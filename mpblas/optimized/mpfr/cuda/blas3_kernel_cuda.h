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
 * MPFR Level 3 BLAS on CUDA: Rsyrk, Rsyr2k, Rsymm, Rtrmm, Rtrsm.
 *
 * Each routine is written once against an arithmetic engine E (fixed
 * precision: cu_fp::cu_freal<PB>, below; runtime precision: cu_mpfr,
 * blas3_rt_cuda.cu) and performs, for every element of the result, the
 * operations of the reference BLAS (mpblas/reference/R*.cpp, which the MPFR
 * libraries use for these routines) in the same order, each correctly
 * rounded to the working precision.  The work is split into independent
 * "lanes":
 *
 *   Rsyrk, Rsyr2k   one element of the referenced triangle of C
 *   Rsymm, Rtrmm, Rtrsm with side = L
 *                   one column of the result (the column's operations are
 *                   those of the reference loops for that column)
 *   Rsymm, Rtrmm, Rtrsm with side = R
 *                   one row of the result
 *
 * Engine interface (all __host__ __device__):
 *   typedef ... V;  struct M { ... };          value / matrix handle
 *   V reg(int slot);                            a register (slot < 8)
 *   V at(const M &, long i, long j);            read-only value
 *   void put(M &, long i, long j, const V &);
 *   void set(V &, const V &); void zero(V &); void one(V &);
 *   void add(V &, const V &, const V &);  sub, mul, div
 *   bool iszero(const V &); bool isone(const V &);
 *   static const bool has_div;
 *
 * Operands are packed column-major with leading dimension = rows.  Each
 * statement of the reference is decomposed into its correctly rounded
 * operations: "c += a*t1 + b*t2" is c = c + ((a*t1) + (b*t2)) and
 * "x = p*q + r*s + u*v" is x = ((p*q) + (r*s)) + (u*v), as gmpfrxx_mkII
 * evaluates them (built without GMPFRXX_MKII_ENABLE_FMA).
 *
 * This header must not include gmpfrxx_mkII or MPFR: it is compiled by nvcc.
 */

#ifndef _MPLAPACK_MPFR_BLAS3_KERNEL_CUDA_H_
#define _MPLAPACK_MPFR_BLAS3_KERNEL_CUDA_H_

#include "Rgemm_kernel_cuda.h"

#if defined(__CUDACC__) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
#define L3_HD __host__ __device__
#else
#define L3_HD
#endif

namespace mplapack_mpfr_cuda {

enum { L3_SYRK = 1, L3_SYR2K = 2, L3_SYMM = 3, L3_TRMM = 4, L3_TRSM = 5 };

struct l3_shape {
    int op;
    int left;     // side == "L" (symm, trmm, trsm)
    int upper;    // uplo == "U"
    int trans;    // trans / transa != "N" (syrk, syr2k, trmm, trsm)
    int nounit;   // diag == "N" (trmm, trsm)
    long m, n, k; // syrk/syr2k: n, k; others: m, n (the size of the result)
    long lda, ldb, ldc;
    int alpha_one, beta_zero, beta_one;
};

// Number of lanes of a call.
L3_HD inline long l3_lanes(const l3_shape &s)
{
    if (s.op == L3_SYRK || s.op == L3_SYR2K)
        return s.n * s.n;
    return s.left ? s.n : s.m;
}

template <class E> struct l3_mats {
    typename E::M al, be, a, b, c; // alpha and beta are 1 x 1
};

// ---- Rsyrk / Rsyr2k: element (i, j) of the referenced triangle ----

template <class E> L3_HD void l3_syrk(E &e, const l3_shape &s, l3_mats<E> &x, long i, long j)
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

template <class E> L3_HD void l3_syr2k(E &e, const l3_shape &s, l3_mats<E> &x, long i, long j)
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

template <class E> L3_HD void l3_symm(E &e, const l3_shape &s, l3_mats<E> &x, long lane)
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
template <class E> L3_HD void l3_axpy1(E &e, typename E::M &B, long i, long j, const typename E::V &t, const typename E::V &w, typename E::V &y, typename E::V &p)
{
    e.set(y, e.at(B, i, j));
    e.mul(p, t, w);
    e.add(y, y, p);
    e.put(B, i, j, y);
}

// b(i,j) = b(i,j) - t * w
template <class E> L3_HD void l3_axmy1(E &e, typename E::M &B, long i, long j, const typename E::V &t, const typename E::V &w, typename E::V &y, typename E::V &p)
{
    e.set(y, e.at(B, i, j));
    e.mul(p, t, w);
    e.sub(y, y, p);
    e.put(B, i, j, y);
}

// b(i,j) = t * b(i,j)
template <class E> L3_HD void l3_scale1(E &e, typename E::M &B, long i, long j, const typename E::V &t, typename E::V &y)
{
    e.set(y, e.at(B, i, j));
    e.mul(y, t, y);
    e.put(B, i, j, y);
}

// ---- Rtrmm: B := alpha*op(A)*B (side L, column lane) or alpha*B*op(A) (side R, row lane) ----
// B is x.c.

template <class E> L3_HD void l3_trmm(E &e, const l3_shape &s, l3_mats<E> &x, long lane)
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

template <class E> L3_HD void l3_trsm(E &e, const l3_shape &s, l3_mats<E> &x, long lane)
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

template <class E> L3_HD void l3_lane(E &e, const l3_shape &s, l3_mats<E> &x, long lane)
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
    default:
        break;
    }
}

// ---- fixed-precision engine (cu_freal<PB>) ----
// No correctly rounded division: Rtrsm with diag = N is not computed here.

template <int PB> struct l3_fixed_engine {
    typedef real_t<PB> V;
    struct M {
        V *p;
        long ld;
    };
    static const bool has_div = false;
    L3_HD V reg(int) const
    {
        V v;
        v.set_zero();
        return v;
    }
    L3_HD V at(const M &a, long i, long j) const { return a.p[i + j * a.ld]; }
    L3_HD void put(M &a, long i, long j, const V &v) const { a.p[i + j * a.ld] = v; }
    L3_HD void set(V &r, const V &a) const { r = a; }
    L3_HD void zero(V &r) const { r.set_zero(); }
    L3_HD void one(V &r) const { r = V::from_double(1.0); }
    L3_HD void add(V &r, const V &a, const V &b) const { r = a + b; }
    L3_HD void sub(V &r, const V &a, const V &b) const { r = a - b; }
    L3_HD void mul(V &r, const V &a, const V &b) const { r = a * b; }
    L3_HD void div(V &r, const V &, const V &) const { r.set_zero(); } // never called
    L3_HD bool iszero(const V &a) const { return a.is_zero(); }
    L3_HD bool isone(const V &a) const { return cu_fp::cu_cmp<PB>(a, V::from_double(1.0)) == 0; }
};

// Level 3 call on packed host operands with the fixed-precision engine.
// alpha/beta are values; A, B are read, C (Rtrmm/Rtrsm: B) is updated.
// Returns 0 on success.  Defined in blas3_device_cuda.cu.
template <int PB> int blas3_fixed(const l3_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, size_t countA, const real_t<PB> *B, size_t countB, real_t<PB> *C, size_t countC);

} // namespace mplapack_mpfr_cuda

#endif
