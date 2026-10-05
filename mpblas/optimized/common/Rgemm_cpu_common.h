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
 * CPU Rgemm driver shared by the dd, qd and binary128 optimized libraries:
 * chooses between the conventional blocked engine (bit-identical to
 * openmp/Rgemm_*_omp.cpp), the Winograd variant of Strassen's algorithm and
 * the Ozaki scheme, from environment variables read at every call
 * (<P> is MPLAPACK_DD, MPLAPACK_QD or MPLAPACK_BINARY128):
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
 */

#ifndef _MPLAPACK_RGEMM_CPU_COMMON_H_
#define _MPLAPACK_RGEMM_CPU_COMMON_H_

#include <string>
#include <vector>
#include "Rgemm_blocked_common.h"
#include "Rgemm_winograd_common.h"
#include "Rgemm_ozaki_common.h"

namespace mplapack_gemm {

inline long env_param(const char *prefix, const char *name, long dflt) { return env_long((std::string(prefix) + name).c_str(), dflt); }

// C = A*B with the blocked engine (Winograd leaves).
template <class K> struct blocked_leaf {
    typedef typename K::value_type T;
    void product(long m, long n, long k, T *C, long ldc, const T *A, long lda, const T *B, long ldb) {
        if ((double)m * n * k >= 32768.0) {
            shape s = {true, true, m, n, k, lda, ldb, ldc};
            gemm_blocked<K>(s, T(1.0), T(0.0), A, B, C);
            return;
        }
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++) {
                T c = 0.0;
                for (long l = 0; l < k; l++)
                    c += B[l + j * ldb] * A[i + l * lda];
                C[i + j * ldc] = c;
            }
    }
};

template <class K> bool gemm_winograd_cpu(const shape &s, const typename K::value_type &alpha, const typename K::value_type &beta, const typename K::value_type *A, const typename K::value_type *B, typename K::value_type *C, long cutoff) {
    typedef typename K::value_type T;
    const long m = s.m, n = s.n, k = s.k;
    std::vector<T> Ab, Bb, P((size_t)m * n);
    const T *Ap = A, *Bp = B;
    long lda = s.lda, ldb = s.ldb;
    if (!s.nota) {
        Ab.resize((size_t)m * k);
        for (long l = 0; l < k; l++)
            for (long i = 0; i < m; i++)
                Ab[i + l * m] = A[l + i * s.lda];
        Ap = Ab.data();
        lda = m;
    }
    if (!s.notb) {
        Bb.resize((size_t)k * n);
        for (long j = 0; j < n; j++)
            for (long l = 0; l < k; l++)
                Bb[l + j * k] = B[j + l * s.ldb];
        Bp = Bb.data();
        ldb = k;
    }
    blocked_leaf<K> leaf;
    winograd_cpu_backend<T, blocked_leaf<K> > be(leaf);
    winograd<T>(be, m, n, k, Ap, lda, Bp, ldb, P.data(), m, cutoff);
    if (be.failed())
        return false;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++) {
            T &c = C[i + j * s.ldc];
            T r = alpha * P[i + j * m];
            if (beta != 0.0)
                r += (beta == 1.0) ? c : beta * c;
            c = r;
        }
    return true;
}

// Returns false when the call is left to the caller (small sizes or the
// blocked engine disabled).  alpha != 0 and m, n > 0 are assumed.
template <class K> bool gemm_cpu(const char *prefix, const shape &s, const typename K::value_type &alpha, const typename K::value_type &beta, const typename K::value_type *A, const typename K::value_type *B, typename K::value_type *C) {
    typedef typename K::value_type T;
    if (s.k > 0 && env_param(prefix, "_GEMM_OZAKI", 0) > 0) {
        if (gemm_ozaki_cpu<T>(s, alpha, beta, A, B, C, (int)env_param(prefix, "_GEMM_OZAKI_SLICES", 0)))
            return true;
    }
    const long cutoff = env_param(prefix, "_GEMM_WINOGRAD_CUTOFF", 0);
    if (cutoff > 0 && s.m > cutoff && s.n > cutoff && s.k > cutoff) {
        if (gemm_winograd_cpu<K>(s, alpha, beta, A, B, C, cutoff))
            return true;
    }
    if (env_param(prefix, "_GEMM_BLOCKED", 1) == 0)
        return false;
    if ((double)s.m * s.n * s.k < (double)env_param(prefix, "_GEMM_MIN_MNK", 32768))
        return false;
    gemm_blocked<K>(s, alpha, beta, A, B, C);
    return true;
}

} // namespace mplapack_gemm

#endif
