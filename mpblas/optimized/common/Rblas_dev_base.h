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
 * What the device text of the reference BLAS routines (blas_dev/) needs
 * besides its element types: the function qualifier and the functions the
 * reference text calls (Mlsame, Mxerbla, max, min).  The element types are
 * Rblas_dev_types.h (dd/td/qd) and mpfr/cuda/Cblas_mpfr_dev_types.h (MPFR).
 */

#ifndef _MPLAPACK_RBLAS_DEV_BASE_H_
#define _MPLAPACK_RBLAS_DEV_BASE_H_

#if defined(__CUDACC__)
#define MPLAPACK_DEV_HD __host__ __device__ __forceinline__
#elif defined(__GNUC__)
#define MPLAPACK_DEV_HD inline __attribute__((always_inline))
#else
#define MPLAPACK_DEV_HD inline
#endif

namespace mplapack_blas_dev {

MPLAPACK_DEV_HD bool lsame(const char *a, const char *b)
{
    char x = a[0], y = b[0];
    if (x >= 'a' && x <= 'z')
        x = (char)(x - 'a' + 'A');
    if (y >= 'a' && y <= 'z')
        y = (char)(y - 'a' + 'A');
    return x == y;
}
MPLAPACK_DEV_HD void xerbla(const char *, int) {} // the host checks the arguments
template <class T> MPLAPACK_DEV_HD T max(T a, T b) { return a > b ? a : b; }
template <class T> MPLAPACK_DEV_HD T min(T a, T b) { return a < b ? a : b; }

// The element types of a device routine: DR (REAL) and DC (COMPLEX) for the
// tag Ops; defined for the dd/td/qd ops in Rblas_dev_types.h and specialized
// for MPFR in mpfr/cuda/Cblas_mpfr_dev_types.h.
template <class Ops> struct cgpu_types;

} // namespace mplapack_blas_dev

#endif
