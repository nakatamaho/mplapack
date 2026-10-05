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
 * Ozaki scheme for Rgemm (dd, qd and binary128 on the CPU; the splitting is
 * shared with the CUDA dd/qd implementation).
 *
 * Every row i of op(A) is scaled by 2^-e_i and every column j of op(B) by
 * 2^-f_j so that all scaled entries are below 1 in magnitude; the scaled
 * entries are then cut into S signed "slices" of beta bits on a common
 * fixed-point grid:
 *
 *   op(A)(i,l) ~ 2^e_i  sum_{s<S} a_s(i,l) 2^-(s+1)beta,   |a_s| < 2^beta,
 *   op(B)(l,j) ~ 2^f_j  sum_{t<S} b_t(l,j) 2^-(t+1)beta,   |b_t| < 2^beta,
 *
 * with beta = floor((53 - ceil(log2 k)) / 2), so that every product
 * D_st = a_s * b_t (a k-term integer dot product per element) is computed
 * exactly by an ordinary double-precision matrix multiplication, whatever
 * the order of its operations.  The products with s + t < S are added in the
 * target precision:
 *
 *   P(i,j) = 2^(e_i + f_j) sum_{g = S-1 .. 0} sum_{s = 0 .. g} D_{s,g-s}(i,j) 2^-(g+2)beta
 *   C := alpha*P + beta*C.
 *
 * S is chosen so that S*beta >= precision + 8 bits.  The error comes from
 * the slices not taken (truncation at 2^-S*beta relative to the largest
 * entry of the row / column) and from the additions in the target
 * precision; it is bounded by about k S 2^(4 - S*beta) 2^(e_i + f_j) per
 * element, i.e. normwise rather than elementwise: entries much smaller than
 * the largest of their row or column lose relative accuracy.
 *
 * Eligibility: all referenced entries of A and B are finite and every nonzero |op(A)(i,l)|, |op(B)(l,j)| lies in
 * [2^-960, 2^960].
 */

#ifndef _MPLAPACK_RGEMM_OZAKI_COMMON_H_
#define _MPLAPACK_RGEMM_OZAKI_COMMON_H_

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include "Rgemm_blocked_common.h"
#include "Rgemm_simd_common.h"

namespace mplapack_gemm {

// Per-type operations (specialized by the users of this header):
//   static const int precision;      // significand bits of T
//   static const int ncomp;          // doubles needed to hold any T exactly
//   static double to_double(const T &);   // nearest double
//   static void add_double(T &acc, double d);   // acc += d
template <class T> struct ozaki_traits;

// For dd_real / qd_real (x[0] is the leading double).
template <class T, int PREC, int NCOMP> struct qd_like_ozaki_traits {
    static const int precision = PREC;
    static const int ncomp = NCOMP;
    static double to_double(const T &x) { return x.x[0]; }
    static void add_double(T &acc, double d) { acc += d; }
};

struct ozaki_plan {
    long m, n, k;
    int S, beta;
    std::vector<int> ea, fb;   // row / column exponents
    std::vector<int32_t> ad;   // a_s(i,l) at s*m*k + i + l*m   (column major, m x k per slice)
    std::vector<int32_t> bd;   // b_t(l,j) at t*k*n + l + j*k   (column major, k x n per slice)
};

inline int ozaki_beta(long k) {
    int lk = 0;
    while ((1L << lk) < k)
        lk++;
    int b = (53 - lk) / 2;
    return b > 26 ? 26 : b;
}

// x * 2^e for |e| <= 2000, exact unless the result over/underflows.
template <class T> inline T ozaki_scale(const T &x, int e) {
    T r = x;
    while (e > 1000) {
        r = r * T(std::ldexp(1.0, 1000));
        e -= 1000;
    }
    while (e < -1000) {
        r = r * T(std::ldexp(1.0, -1000));
        e += 1000;
    }
    return r * T(std::ldexp(1.0, e));
}

// Signed fixed-point integer of up to 384 bits (two's complement limbs).
struct ozaki_fixed {
    static const int L = 6;
    uint64_t w[L];
    void clear() { std::memset(w, 0, sizeof(w)); }
    // += sign * mant * 2^sh, mant < 2^53; bits below 2^0 are dropped
    // (the magnitude is truncated).
    void add(uint64_t mant, int sh, bool neg) {
        uint64_t t[L];
        std::memset(t, 0, sizeof(t));
        if (sh < 0) {
            if (sh <= -64)
                return;
            mant >>= -sh;
            sh = 0;
        }
        const int q = sh / 64, r = sh % 64;
        if (q < L)
            t[q] = mant << r;
        if (r != 0 && q + 1 < L)
            t[q + 1] = mant >> (64 - r);
        if (neg) {
            // two's complement of t
            uint64_t carry = 1;
            for (int i = 0; i < L; i++) {
                const uint64_t v = ~t[i] + carry;
                carry = (carry && v == 0) ? 1 : 0;
                t[i] = v;
            }
        }
        uint64_t carry = 0;
        for (int i = 0; i < L; i++) {
            const uint64_t a = w[i], s1 = a + t[i], s2 = s1 + carry;
            carry = (s1 < a) || (s2 < s1);
            w[i] = s2;
        }
    }
    bool negative() const { return (w[L - 1] >> 63) != 0; }
    void negate() {
        uint64_t carry = 1;
        for (int i = 0; i < L; i++) {
            const uint64_t v = ~w[i] + carry;
            carry = (carry && v == 0) ? 1 : 0;
            w[i] = v;
        }
    }
    // bits [lo, lo + n) of a nonnegative value, n <= 32
    uint32_t bits(int lo, int n) const {
        const int q = lo / 64, r = lo % 64;
        uint64_t v = w[q] >> r;
        if (r != 0 && q + 1 < L)
            v |= w[q + 1] << (64 - r);
        return (uint32_t)(v & ((1ULL << n) - 1));
    }
};

// Signed digits of y (|y| < 1) on the grid 2^-W, W = S*beta: y ~ sum_s d[s] 2^-(s+1)beta.
template <class T> inline void ozaki_digits(const T &y, int S, int beta, int32_t *d, long stride) {
    typedef ozaki_traits<T> tr;
    const int W = S * beta;
    ozaki_fixed v;
    v.clear();
    T r = y;
    for (int c = 0; c < tr::ncomp; c++) {
        const double x = tr::to_double(r);
        if (x == 0.0)
            break;
        int q;
        const double f = std::frexp(std::fabs(x), &q); // |x| = f 2^q, f in [0.5, 1)
        const uint64_t mant = (uint64_t)std::ldexp(f, 53);
        v.add(mant, q - 53 + W, x < 0.0);
        r = r - T(x);
    }
    const bool neg = v.negative();
    if (neg)
        v.negate();
    for (int s = 0; s < S; s++) {
        const int32_t dv = (int32_t)v.bits(W - (s + 1) * beta, beta);
        d[s * stride] = neg ? -dv : dv;
    }
}

// Exponent e with |x| < 2^e for every x of the set whose largest magnitude
// (as a double) is mx; 0 when mx == 0.
inline int ozaki_exponent(double mx) { return mx == 0.0 ? 0 : std::ilogb(mx) + 2; }

inline bool ozaki_in_range(double x) { return x == 0.0 || (std::fabs(x) >= std::ldexp(1.0, -960) && std::fabs(x) <= std::ldexp(1.0, 960)); }

// Splits op(A) and op(B).  Returns false when the scheme is not applicable.
template <class T> bool ozaki_split(const shape &s, const T *A, const T *B, int S, ozaki_plan &p) {
    typedef ozaki_traits<T> tr;
    const long m = s.m, n = s.n, k = s.k;
    p.m = m;
    p.n = n;
    p.k = k;
    p.beta = ozaki_beta(k);
    p.S = S > 0 ? S : (tr::precision + 8 + p.beta - 1) / p.beta;
    if (p.S * p.beta > 320) // capacity of ozaki_fixed
        p.S = 320 / p.beta;
    p.ea.assign(m, 0);
    p.fb.assign(n, 0);
    bool ok = true;
#ifdef _OPENMP
#pragma omp parallel for reduction(&& : ok)
#endif
    for (long i = 0; i < m; i++) {
        double mx = 0.0;
        for (long l = 0; l < k; l++) {
            const double x = tr::to_double(opA(s, A, i, l));
            if (!std::isfinite(x) || !ozaki_in_range(x))
                ok = false;
            else if (std::fabs(x) > mx)
                mx = std::fabs(x);
        }
        p.ea[i] = ozaki_exponent(mx);
    }
#ifdef _OPENMP
#pragma omp parallel for reduction(&& : ok)
#endif
    for (long j = 0; j < n; j++) {
        double mx = 0.0;
        for (long l = 0; l < k; l++) {
            const double x = tr::to_double(opB(s, B, l, j));
            if (!std::isfinite(x) || !ozaki_in_range(x))
                ok = false;
            else if (std::fabs(x) > mx)
                mx = std::fabs(x);
        }
        p.fb[j] = ozaki_exponent(mx);
    }
    if (!ok)
        return false;
    const int SS = p.S, beta = p.beta;
    p.ad.assign((size_t)SS * m * k, 0);
    p.bd.assign((size_t)SS * k * n, 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long l = 0; l < k; l++)
        for (long i = 0; i < m; i++)
            ozaki_digits(ozaki_scale(opA(s, A, i, l), -p.ea[i]), SS, beta, &p.ad[i + l * m], m * k);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long j = 0; j < n; j++)
        for (long l = 0; l < k; l++)
            ozaki_digits(ozaki_scale(opB(s, B, l, j), -p.fb[j]), SS, beta, &p.bd[l + j * k], k * n);
    return true;
}

// ---- exact double-precision product of integer matrices ----

// C (m x n) = A (m x k) * B (k x n), all column major, A and B int32 with
// |sum| < 2^53 so that every operation is exact.
MPLAPACK_GEMM_CLONES static void ozaki_dgemm_block(long mc, long nc, long kc, const double *ap, const double *bp, double *C, long ldc, bool first) {
    const int MR = 8, NR = 4;
    for (long jr = 0; jr < nc; jr += NR) {
        const long nr = nc - jr < NR ? nc - jr : NR;
        for (long ir = 0; ir < mc; ir += MR) {
            const long mr = mc - ir < MR ? mc - ir : MR;
            double acc[NR][MR];
            for (int jj = 0; jj < NR; jj++)
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] = (!first && jj < nr && ii < mr) ? C[(ir + ii) + (jr + jj) * ldc] : 0.0;
            const double *a = ap + ir * kc;  // [l][MR] for this row panel
            const double *b = bp + jr * kc;  // [l][NR] for this column panel
            for (long l = 0; l < kc; l++) {
                for (int jj = 0; jj < NR; jj++) {
                    const double bv = b[l * NR + jj];
                    MPLAPACK_GEMM_SIMD
                    for (int ii = 0; ii < MR; ii++)
                        acc[jj][ii] += a[l * MR + ii] * bv;
                }
            }
            for (int jj = 0; jj < nr; jj++)
                for (int ii = 0; ii < mr; ii++)
                    C[(ir + ii) + (jr + jj) * ldc] = acc[jj][ii];
        }
    }
}

inline void ozaki_dgemm(long m, long n, long k, const int32_t *A, long lda, const int32_t *B, long ldb, double *C, long ldc) {
    const long MC = 128, NC = 256, KC = 256, MR = 8, NR = 4;
    const long nbi = (m + MC - 1) / MC, nbj = (n + NC - 1) / NC;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<double> ap((size_t)(MC + MR) * KC), bp((size_t)(NC + NR) * KC);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (long blk = 0; blk < nbi * nbj; blk++) {
            const long i0 = (blk % nbi) * MC, j0 = (blk / nbi) * NC;
            const long mc = m - i0 < MC ? m - i0 : MC, nc = n - j0 < NC ? n - j0 : NC;
            if (k == 0) {
                for (long j = 0; j < nc; j++)
                    for (long i = 0; i < mc; i++)
                        C[(i0 + i) + (j0 + j) * ldc] = 0.0;
                continue;
            }
            for (long p0 = 0; p0 < k; p0 += KC) {
                const long kc = k - p0 < KC ? k - p0 : KC;
                // A panel: row panels of MR, [l][MR], zero padded
                for (long ir = 0; ir < mc; ir += MR)
                    for (long l = 0; l < kc; l++)
                        for (long ii = 0; ii < MR; ii++)
                            ap[ir * kc + l * MR + ii] = (ir + ii < mc) ? (double)A[(i0 + ir + ii) + (p0 + l) * lda] : 0.0;
                for (long jr = 0; jr < nc; jr += NR)
                    for (long l = 0; l < kc; l++)
                        for (long jj = 0; jj < NR; jj++)
                            bp[jr * kc + l * NR + jj] = (jr + jj < nc) ? (double)B[(p0 + l) + (j0 + jr + jj) * ldb] : 0.0;
                ozaki_dgemm_block(mc, nc, kc, ap.data(), bp.data(), C + i0 + j0 * ldc, ldc, p0 == 0);
            }
        }
    }
}

// acc(i,j) += D(i,j) 2^-sh for the plan's element order
template <class T> inline void ozaki_accumulate(long m, long n, T *acc, const double *D, int sh) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            ozaki_traits<T>::add_double(acc[i + j * m], std::ldexp(D[i + j * m], -sh));
}

// C(i,j) := alpha * 2^(e_i + f_j) acc(i,j) + beta*C(i,j)
template <class T> inline T ozaki_finish(const T &accv, int e, int f, const T &alpha, const T &beta, const T &c) {
    T p = (accv * T(std::ldexp(1.0, e))) * T(std::ldexp(1.0, f));
    T r = alpha * p;
    if (beta != 0.0)
        r += (beta == 1.0) ? c : beta * c;
    return r;
}

// The whole scheme on the CPU.  Returns false (C untouched) when not applicable.
template <class T> bool gemm_ozaki_cpu(const shape &s, const T &alpha, const T &beta, const T *A, const T *B, T *C, int S) {
    ozaki_plan p;
    if (!ozaki_split(s, A, B, S, p))
        return false;
    const long m = s.m, n = s.n, k = s.k;
    std::vector<T> acc((size_t)m * n, T(0.0));
    std::vector<double> D((size_t)m * n);
    for (int g = p.S - 1; g >= 0; g--)
        for (int a = 0; a <= g; a++) {
            const int b = g - a;
            ozaki_dgemm(m, n, k, &p.ad[(size_t)a * m * k], m, &p.bd[(size_t)b * k * n], k, D.data(), m);
            ozaki_accumulate(m, n, acc.data(), D.data(), (g + 2) * p.beta);
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            C[i + j * s.ldc] = ozaki_finish(acc[i + j * m], p.ea[i], p.fb[j], alpha, beta, C[i + j * s.ldc]);
    return true;
}

} // namespace mplapack_gemm

#endif
