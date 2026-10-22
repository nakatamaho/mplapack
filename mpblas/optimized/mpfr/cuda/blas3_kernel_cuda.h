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
 * MPFR BLAS on CUDA: the lane routines of ../../common/Rblas_lanes_common.h
 * (Level 3: Rsyrk, Rsyr2k, Rsymm, Rtrmm, Rtrsm; Level 2 and 1) run with an
 * MPFR-exact engine, fixed precision cu_fp::cu_freal<PB> (below) or runtime
 * precision cu_mpfr (Rgemm_rt_cuda.cu).  Every operation is correctly
 * rounded to the working precision, as in the reference BLAS, so the result
 * is the reference one bit for bit (cu_freal: zeros are +0).
 *
 * Operands are packed column-major with leading dimension = rows (vectors:
 * increment 1).  This header must not include gmpfrxx_mkII or MPFR: it is
 * compiled by nvcc.
 */

#ifndef _MPLAPACK_MPFR_BLAS3_KERNEL_CUDA_H_
#define _MPLAPACK_MPFR_BLAS3_KERNEL_CUDA_H_

#include "Rgemm_kernel_cuda.h"
#include "../../common/Rblas_lanes_common.h"

#define L3_HD MPLAPACK_BLAS_HD

namespace mplapack_mpfr_cuda {

using namespace mplapack_blas;
typedef blas_shape l3_shape;

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
template <int PB> int blas3_fixed(const blas_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, size_t countA, const real_t<PB> *B, size_t countB, real_t<PB> *C, size_t countC);

} // namespace mplapack_mpfr_cuda

#endif
