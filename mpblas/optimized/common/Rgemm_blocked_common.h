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
 * The shape of an Rgemm call and the helpers shared by the Winograd, Ozaki
 * and CUDA Rgemm code.  The blocked CPU engine is Rblas_block_common.h.
 */

#ifndef _MPLAPACK_RGEMM_BLOCKED_COMMON_H_
#define _MPLAPACK_RGEMM_BLOCKED_COMMON_H_

#include <cstdlib>
#include <cstring>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mplapack_gemm {

// C (m x n) := alpha*op(A)*op(B) + beta*C, column major.
struct shape {
    bool nota, notb;
    long m, n, k, lda, ldb, ldc;
};

// Integer value of an environment variable, or dflt when unset or empty.
inline long env_long(const char *name, long dflt) {
    const char *v = std::getenv(name);
    if (v == NULL || *v == '\0')
        return dflt;
    return std::strtol(v, NULL, 10);
}

template <class T> inline const T &opA(const shape &s, const T *A, long i, long l) { return s.nota ? A[i + l * s.lda] : A[l + i * s.lda]; }
template <class T> inline const T &opB(const shape &s, const T *B, long l, long j) { return s.notb ? B[l + j * s.ldb] : B[j + l * s.ldb]; }

// The value C(i,j) takes before the products are added (the first loop of
// Rgemm_*_omp).
template <class T> inline T scaled_c(const T &beta, const T &c) {
    if (beta == 0.0)
        return T(0.0);
    if (beta != 1.0)
        return beta * c;
    return c;
}

} // namespace mplapack_gemm

#endif
