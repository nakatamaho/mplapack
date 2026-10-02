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
 * Host implementation of mplapack_mpfr_cuda::gemm<PB>.  It runs the same
 * element routines as the CUDA kernels (Rgemm_kernel_cuda.h) on the CPU, so the
 * GPU code path can be validated against MPFR on machines without a GPU.
 * It is linked instead of Rgemm_device_cuda.cu by the host-emulation test.
 */

#include <vector>
#include "Rgemm_kernel_cuda.h"

namespace mplapack_mpfr_cuda {

template <int PB> int gemm(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C)
{
    std::vector<real_t<PB>> AB(s.nota ? (size_t)s.k * s.n : 0);
    if (s.nota) {
        for (long j = 0; j < s.n; j++)
            for (long l = 0; l < s.k; l++)
                gemm_alpha_opB<PB>(s, alpha, B, AB.data(), l, j);
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long j = 0; j < s.n; j++)
        for (long i = 0; i < s.m; i++)
            gemm_element<PB>(s, alpha, beta, A, B, AB.data(), C, i, j);
    return 0;
}

int device_count() { return 1; }

template int gemm<512>(const gemm_shape &, const real_t<512> &, const real_t<512> &, const real_t<512> *, const real_t<512> *, real_t<512> *);
template int gemm<1024>(const gemm_shape &, const real_t<1024> &, const real_t<1024> &, const real_t<1024> *, const real_t<1024> *, real_t<1024> *);

} // namespace mplapack_mpfr_cuda
