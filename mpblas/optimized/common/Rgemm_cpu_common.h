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
 * CPU Rgemm driver shared by the dd, td, qd and binary128 optimized libraries:
 * chooses between the conventional blocked engine (bit-identical to
 * openmp/Rgemm_*_omp.cpp), the Winograd variant of Strassen's algorithm and
 * the Ozaki scheme, from environment variables read at every call
 * (<P> is MPLAPACK_DD, MPLAPACK_TD, MPLAPACK_QD or MPLAPACK_BINARY128):
 *
 *   <P>_GEMM_BLOCKED=0            use the plain OpenMP loops (Rgemm_*_omp)
 *   <P>_GEMM_MIN_MNK=N            smallest m*n*k for the blocked engine (32768)
 *   <P>_GEMM_WINOGRAD_CUTOFF=N    Winograd when min(m,n,k) > N (0: off)
 *   <P>_GEMM_OZAKI=1              Ozaki scheme (off)
 *   <P>_GEMM_OZAKI_SLICES=S       number of slices (default: precision+8 bits)
 *
 * Winograd and Ozaki are opt-in: their results are not those of the
 * conventional Rgemm and their errors are bounded normwise only.  Ozaki takes
 * precedence over Winograd and falls back to the other paths when it is not
 * applicable (non-finite or out-of-range entries).
 *
 * The other libraries (double, binary80, MPFR, GMP; <P> = MPLAPACK_DOUBLE,
 * MPLAPACK_BINARY80, MPLAPACK_MPFR, MPLAPACK_GMP) take only the Winograd
 * variant (gemm_winograd_env, MPLAPACK_GEMM_NO_OZAKI before the include),
 * with the leaves of their own Rgemm, and so do the complex Cgemm of every
 * library (<P>_CGEMM_WINOGRAD_CUTOFF).  The temporaries are made and
 * computed in threads with the caller's context (the MPFR and GMP default
 * precision), the packed copies keep the precision of their sources.
 */

#ifndef _MPLAPACK_RGEMM_CPU_COMMON_H_
#define _MPLAPACK_RGEMM_CPU_COMMON_H_

#include <string>
#include <vector>
#include "Rgemm_blocked_common.h"
#include "Rblas_cpu_common.h"
#include "Rgemm_winograd_common.h"
#ifndef MPLAPACK_GEMM_NO_OZAKI
#include "Rgemm_ozaki_common.h"
#endif

namespace mplapack_gemm {

inline long env_param(const char *prefix, const char *name, long dflt) { return env_long((std::string(prefix) + name).c_str(), dflt); }

// C = A*B with the small-block engine (Winograd leaves); where the engine is
// not used for T (binary128), Rgemm of the library (the leaves are at most
// the cutoff in one dimension, so it does not recurse).
template <class T> struct blocked_leaf {
    void product(long m, long n, long k, T *C, long ldc, const T *A, long lda, const T *B, long ldb) {
        if (mplapack_blas::blk_gemm<mplapack_blas::blas_thread_context>(true, true, m, n, k, T(1.0), A, lda, B, ldb, T(0.0), C, ldc, 0.0))
            return;
        const T one = 1.0, zero = 0.0;
        Rgemm("N", "N", (INTEGER)m, (INTEGER)n, (INTEGER)k, one, const_cast<T *>(A), (INTEGER)lda, const_cast<T *>(B), (INTEGER)ldb, zero, C, (INTEGER)ldc);
    }
};

// conj(x) for the complex types, x otherwise.
template <class T> inline T winograd_conj(const T &x) {
    if constexpr (mplapack_blas::blk_is_cplx<T>::value)
        return T(conj(x));
    else
        return x;
}

// The Winograd backend of the CPU: element loops over the columns in
// threads with the caller's context; zero is the zero of F (complex types
// pass their own).
template <class F, class Leaf> class winograd_ctx_backend {
  public:
    winograd_ctx_backend(Leaf &leaf, const F &zero) : leaf_(leaf), zero_(zero), failed_(false) {}
    F *alloc(size_t count) {
        F *p = new (std::nothrow) F[count ? count : 1];
        if (p == NULL)
            failed_ = true;
        return p;
    }
    void release(F *p) { delete[] p; }
    bool failed() const { return failed_; }
    void zero(long r, long c, F *X, long ldx) {
        cols(r, c, [&](long j) {
            for (long i = 0; i < r; i++)
                X[i + j * ldx] = zero_;
        });
    }
    void copy(long r, long c, F *D, long ldd, const F *S, long lds) {
        cols(r, c, [&](long j) {
            for (long i = 0; i < r; i++)
                D[i + j * ldd] = S[i + j * lds];
        });
    }
    void addsub(long r, long c, F *Z, long ldz, const F *X, long ldx, const F *Y, long ldy, int sub) {
        cols(r, c, [&](long j) {
            for (long i = 0; i < r; i++)
                if (sub)
                    Z[i + j * ldz] = X[i + j * ldx] - Y[i + j * ldy];
                else
                    Z[i + j * ldz] = X[i + j * ldx] + Y[i + j * ldy];
        });
    }
    void leaf(long m, long n, long k, F *C, long ldc, const F *A, long lda, const F *B, long ldb) { leaf_.product(m, n, k, C, ldc, A, lda, B, ldb); }

  private:
    template <class G> void cols(long r, long c, G g) { mplapack_blas::blk_columns<mplapack_blas::blas_thread_context>(c, mplapack_blas::blk_parallel((double)r * c * 16), g); }
    Leaf &leaf_;
    F zero_;
    bool failed_;
};

// C := alpha*P + beta*C with P = op(A)*op(B) by the Winograd recursion (leaf:
// C = A*B of the library, conj: op(x) = x^H for a transposed x); op(A) and
// op(B) are packed into copies first.
template <class T, class Leaf> bool gemm_winograd_cpu(const shape &s, const T &alpha, const T &beta, const T *A, const T *B, T *C, long cutoff, Leaf &leaf, const T &zero, const T &one, bool conja = false, bool conjb = false) {
    const long m = s.m, n = s.n, k = s.k;
    std::vector<T> Ab, Bb, P((size_t)m * n);
    const T *Ap = A, *Bp = B;
    long lda = s.lda, ldb = s.ldb;
    if (!s.nota) {
        Ab.reserve((size_t)m * k);
        for (long l = 0; l < k; l++)
            for (long i = 0; i < m; i++)
                Ab.push_back(conja ? winograd_conj(A[l + i * s.lda]) : T(A[l + i * s.lda])); // copies keep the precision
        Ap = Ab.data();
        lda = m;
    }
    if (!s.notb) {
        Bb.reserve((size_t)k * n);
        for (long j = 0; j < n; j++)
            for (long l = 0; l < k; l++)
                Bb.push_back(conjb ? winograd_conj(B[j + l * s.ldb]) : T(B[j + l * s.ldb]));
        Bp = Bb.data();
        ldb = k;
    }
    winograd_ctx_backend<T, Leaf> be(leaf, zero);
    winograd<T>(be, m, n, k, Ap, lda, Bp, ldb, P.data(), m, cutoff);
    if (be.failed())
        return false;
    const bool beta_zero = beta == zero, beta_one = beta == one;
    mplapack_blas::blk_columns<mplapack_blas::blas_thread_context>(n, mplapack_blas::blk_parallel((double)m * n * 64), [&](long j) {
        for (long i = 0; i < m; i++) {
            T &c = C[i + j * s.ldc];
            T r = alpha * P[i + j * m];
            if (beta_one)
                r += c;
            else if (!beta_zero)
                r += beta * c;
            c = r;
        }
    });
    return true;
}

// The Winograd variant when <prefix><var> = N > 0 and min(m,n,k) > N.
template <class T, class Leaf> bool gemm_winograd_env(const char *prefix, const char *var, const shape &s, const T &alpha, const T &beta, const T *A, const T *B, T *C, Leaf &leaf, const T &zero, const T &one, bool conja = false, bool conjb = false) {
    const long cutoff = env_param(prefix, var, 0);
    if (s.k <= 0 || cutoff <= 0 || s.m <= cutoff || s.n <= cutoff || s.k <= cutoff)
        return false;
    return gemm_winograd_cpu<T, Leaf>(s, alpha, beta, A, B, C, cutoff, leaf, zero, one, conja, conjb);
}

#ifndef MPLAPACK_GEMM_NO_OZAKI
// Returns false when the call is left to the caller (small sizes or the
// blocked engine disabled).  alpha != 0 and m, n > 0 are assumed.
template <class T> bool gemm_cpu(const char *prefix, const shape &s, const T &alpha, const T &beta, const T *A, const T *B, T *C) {
    if (s.k > 0 && env_param(prefix, "_GEMM_OZAKI", 0) > 0) {
        if (gemm_ozaki_cpu<T>(s, alpha, beta, A, B, C, (int)env_param(prefix, "_GEMM_OZAKI_SLICES", 0)))
            return true;
    }
    blocked_leaf<T> leaf;
    if (gemm_winograd_env<T>(prefix, "_GEMM_WINOGRAD_CUTOFF", s, alpha, beta, A, B, C, leaf, T(0.0), T(1.0)))
        return true;
    if (env_param(prefix, "_GEMM_BLOCKED", 1) == 0)
        return false;
    if ((double)s.m * s.n * s.k < (double)env_param(prefix, "_GEMM_MIN_MNK", 32768))
        return false;
    return mplapack_blas::blk_gemm<mplapack_blas::blas_thread_context>(s.nota, s.notb, s.m, s.n, s.k, alpha, A, s.lda, B, s.ldb, beta, C, s.ldc, 0.0);
}
#endif

} // namespace mplapack_gemm

#endif
