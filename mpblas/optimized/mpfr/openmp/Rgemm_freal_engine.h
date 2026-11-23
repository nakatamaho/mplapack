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

// The fixed-precision MPFR Rgemm on the CPU (included by Rgemm_freal_x86.cpp
// and Rgemm_freal_generic.cpp, each with its own MPFREAL_NS): the elements
// are copied into cu_freal<PB> (../lgpl/Rfreal_host.h, LGPL-3.0-or-later) tile by tile and
// every operation of openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp is done in the same
// order, rounded to nearest-even -- the result is theirs bit for bit,
// signed zeros included.  Rgemm_freal_omp.cpp has already checked that every
// operand is a regular number or zero of precision PB and that no
// intermediate can leave the exponent range.
//
// The includer has included <cstring>, <mpfr.h> and <mpblas_mpfr.h> and
// defines MPFREAL_NS and MPFREAL_ENTRY (the name of the entry function; the
// fixed-precision LAPACK of Rfreal_blas_inst.h uses freal_gemm on fr<PB>
// arrays and defines none).

#include "Rgemm_freal.h"
#include "../lgpl/Rfreal_host.h"

namespace MPFREAL_NS {

template <int PB> struct freal_gemm {
    typedef cu_freal<PB> F;
    static const int N = F::N;
    // Tiles of C of MB x NB elements, l in panels of KC steps (the packed
    // panels of A and B stay in the L2 cache).
    static const long MB = 32, NB = 32;
    static const long KC = (256L * 1024) / (MB * (long)sizeof(F)) < 16 ? 16 : (256L * 1024) / (MB * (long)sizeof(F));

    static inline void load(F &f, mpfr_srcptr x)
    {
        if (mpfr_zero_p(x)) {
            f.set_zero();
            f.sign = mpfr_signbit(x) ? -1 : 1;
            return;
        }
        f.sign = mpfr_sgn(x) < 0 ? -1 : 1;
        f.exp = mpfr_get_exp(x);
        std::memcpy(f.m, x->_mpfr_d, sizeof f.m);
    }
    static inline void store(mpfr_ptr x, const F &f)
    {
        if (f.is_zero()) {
            mpfr_set_zero(x, f.sign);
            return;
        }
        std::memcpy(x->_mpfr_d, f.m, sizeof f.m);
        x->_mpfr_exp = f.exp;
        x->_mpfr_sign = f.sign;
    }
    // the elements of fr<PB> (Rfreal_types.h: a cu_freal with the sign of its zero)
    static inline void load(F &f, const mpfr_class &x) { load(f, x.get_mpfr_t()); }
    static inline void store(mpfr_class &x, const F &f) { store(x.get_mpfr_t(), f); }
    template <class E> static inline void load(F &f, const E &x) { f = x.v; }
    template <class E> static inline void store(E &x, const F &f) { x.v = f; }
    // MPFR's signed zeros: a product with a zero factor is zero with the
    // product of the signs; a sum of two zeros is -0 only if both are -0.
    static inline F mul(const F &a, const F &b)
    {
        if (a.is_zero() || b.is_zero()) {
            F r;
            r.set_zero();
            r.sign = a.sign * b.sign;
            return r;
        }
        return cu_fmul<PB>(a, b);
    }
    static inline F add(const F &a, const F &b)
    {
        if (a.is_zero() && b.is_zero()) {
            F r;
            r.set_zero();
            r.sign = (a.sign < 0 && b.sign < 0) ? -1 : 1;
            return r;
        }
        return cu_fadd<PB>(a, b);
    }

    // beta C as the OpenMP loops form it: 0 (beta == 0), C (beta == 1), beta*C
    template <class E> static inline F scaled_c(const F &beta, int beta_kind, const E &c)
    {
        F v;
        if (beta_kind == 0) {
            v.set_zero();
            return v;
        }
        load(v, c);
        return beta_kind == 1 ? v : mul(beta, v);
    }

    // E: mpfr_class or fr<PB>
    template <class E>
    static void run(bool nota, bool notb, long m, long n, long k, const E &alpha_, const E *A, long lda, const E *B, long ldb, const E &beta_, int beta_kind, E *C,
                    long ldc, bool parallel)
    {
        F alpha, beta;
        load(alpha, alpha_);
        load(beta, beta_);
        const long mt = (m + MB - 1) / MB, nt = (n + NB - 1) / NB;
#ifdef _OPENMP
#pragma omp parallel if (parallel)
#endif
        {
            F *ap = new F[MB * KC], *bp = new F[KC * NB], *ct = new F[MB * NB];
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(dynamic)
#endif
            for (long jt = 0; jt < nt; jt++)
                for (long it = 0; it < mt; it++) {
                    const long i0 = it * MB, j0 = jt * NB, mb = m - i0 < MB ? m - i0 : MB, nb = n - j0 < NB ? n - j0 : NB;
                    // C := beta C first (NN, NT) / the dot products from 0 (TN, TT)
                    for (long jj = 0; jj < nb; jj++)
                        for (long ii = 0; ii < mb; ii++) {
                            if (nota)
                                ct[ii + jj * MB] = scaled_c(beta, beta_kind, C[(i0 + ii) + (j0 + jj) * ldc]);
                            else
                                ct[ii + jj * MB].set_zero();
                        }
                    for (long l0 = 0; l0 < k; l0 += KC) {
                        const long kc = k - l0 < KC ? k - l0 : KC;
                        // bp[ll + jj KC] = op(B)(l0 + ll, j0 + jj)
                        for (long jj = 0; jj < nb; jj++)
                            for (long ll = 0; ll < kc; ll++)
                                load(bp[ll + jj * KC], (notb ? B[(l0 + ll) + (j0 + jj) * ldb] : B[(j0 + jj) + (l0 + ll) * ldb]));
                        if (nota) {
                            // ap[ii + ll MB] = A(i0 + ii, l0 + ll);
                            // C(i, j) += (alpha op(B)(l, j)) A(i, l), l ascending
                            for (long ll = 0; ll < kc; ll++)
                                for (long ii = 0; ii < mb; ii++)
                                    load(ap[ii + ll * MB], A[(i0 + ii) + (l0 + ll) * lda]);
                            for (long jj = 0; jj < nb; jj++) {
                                F *c = ct + jj * MB;
                                for (long ll = 0; ll < kc; ll++) {
                                    const F temp = mul(alpha, bp[ll + jj * KC]);
                                    const F *a = ap + ll * MB;
                                    for (long ii = 0; ii < mb; ii++)
                                        c[ii] = add(c[ii], mul(temp, a[ii]));
                                }
                            }
                        } else {
                            // ap[ll + ii KC] = A(l0 + ll, i0 + ii);
                            // temp(i, j) += A(l, i) op(B)(l, j), l ascending
                            for (long ii = 0; ii < mb; ii++)
                                for (long ll = 0; ll < kc; ll++)
                                    load(ap[ll + ii * KC], A[(l0 + ll) + (i0 + ii) * lda]);
                            for (long jj = 0; jj < nb; jj++)
                                for (long ii = 0; ii < mb; ii++) {
                                    F t = ct[ii + jj * MB];
                                    const F *a = ap + ii * KC, *b = bp + jj * KC;
                                    for (long ll = 0; ll < kc; ll++)
                                        t = add(t, mul(a[ll], b[ll]));
                                    ct[ii + jj * MB] = t;
                                }
                        }
                    }
                    for (long jj = 0; jj < nb; jj++)
                        for (long ii = 0; ii < mb; ii++) {
                            E &c = C[(i0 + ii) + (j0 + jj) * ldc];
                            if (nota)
                                store(c, ct[ii + jj * MB]);
                            else // C := beta C, then C += alpha temp
                                store(c, add(scaled_c(beta, beta_kind, c), mul(alpha, ct[ii + jj * MB])));
                        }
                }
            delete[] ap;
            delete[] bp;
            delete[] ct;
        }
    }
};

} // namespace MPFREAL_NS

#ifdef MPFREAL_ENTRY
bool MPFREAL_ENTRY(int prec, bool nota, bool notb, long m, long n, long k, const mpfr_class &alpha, const mpfr_class *A, long lda, const mpfr_class *B, long ldb,
                   const mpfr_class &beta, int beta_kind, mpfr_class *C, long ldc, bool parallel)
{
    switch (prec) {
#define MPFREAL_CASE(P)                                                                                                                                        \
    case P:                                                                                                                                                    \
        MPFREAL_NS::freal_gemm<P>::run(nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, beta_kind, C, ldc, parallel);                                         \
        return true;
        MPFREAL_PRECISIONS(MPFREAL_CASE)
#undef MPFREAL_CASE
    default:
        return false;
    }
}
#endif
