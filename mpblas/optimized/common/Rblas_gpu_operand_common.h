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
 * An operand of a GPU BLAS call as the host sees it: the span
 * host[lo .. lo + count - 1] of an operand of leading dimension ld (vectors:
 * the increment); host points at its logical element (0, 0).  download: the
 * GPU writes it.  rows x cols and region (REGION_*) describe the elements
 * the call references: only those are copied to the device and back, so a
 * call never reads or writes the other elements of the span (rows outside a
 * submatrix, the elements between those of a strided vector, the other
 * triangle), which concurrent calls on other parts of the same array may be
 * updating.
 */

#ifndef _MPLAPACK_RBLAS_GPU_OPERAND_COMMON_H_
#define _MPLAPACK_RBLAS_GPU_OPERAND_COMMON_H_

#include <cstring>

namespace mplapack_blas {

enum { REGION_ALL = 0, REGION_UPPER = 1, REGION_LOWER = 2, REGION_NODIAG = 4 };

template <class F> struct gpu_operand {
    const F *host;
    long ld, lo, count;
    bool download;
    long rows, cols;
    int region;
};

template <class F> inline gpu_operand<F> gpu_scalar_operand(const F *x)
{
    gpu_operand<F> o = {x, 1, 0, 1, false, 1, 1, REGION_ALL};
    return o;
}

// The referenced elements of column j: rows [i0, i1).
template <class F> inline void gpu_operand_rows(const gpu_operand<F> &o, long j, long &i0, long &i1)
{
    i0 = 0, i1 = o.rows;
    switch (o.region & 3) {
    case REGION_UPPER:
        i1 = (j + 1 < o.rows) ? j + 1 : o.rows;
        break;
    case REGION_LOWER:
        i0 = j;
        break;
    default:
        break;
    }
    if ((o.region & REGION_NODIAG) && j >= i0 && j < i1) {
        if ((o.region & 3) == REGION_UPPER)
            i1 = j;
        else if ((o.region & 3) == REGION_LOWER)
            i0 = j + 1;
    }
    if (i1 < i0)
        i1 = i0;
}

// True when the referenced elements are the whole span (one block copy).
template <class F> inline bool gpu_operand_dense(const gpu_operand<F> &o)
{
    return o.region == REGION_ALL && (o.cols <= 1 || o.ld == o.rows);
}

// Copy the referenced elements between the host operand and a buffer laid
// out as its span (buf[t] is host[lo + t]).
template <class F> inline void gpu_operand_gather(const gpu_operand<F> &o, F *buf)
{
    for (long j = 0; j < o.cols; j++) {
        long i0, i1;
        gpu_operand_rows(o, j, i0, i1);
        if (i1 > i0)
            std::memcpy(buf + (j * o.ld + i0 - o.lo), o.host + (j * o.ld + i0), (size_t)(i1 - i0) * sizeof(F));
    }
}

template <class F> inline void gpu_operand_scatter(const gpu_operand<F> &o, const F *buf)
{
    F *h = const_cast<F *>(o.host);
    for (long j = 0; j < o.cols; j++) {
        long i0, i1;
        gpu_operand_rows(o, j, i0, i1);
        if (i1 > i0)
            std::memcpy(h + (j * o.ld + i0), buf + (j * o.ld + i0 - o.lo), (size_t)(i1 - i0) * sizeof(F));
    }
}

} // namespace mplapack_blas

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#include <vector>

namespace mplapack_blas {

// dev holds the span of o (dev[t] is host[lo + t]).  Rectangles with ld > 0
// use one 2D copy; triangles and negative increments go through a host
// staging buffer.
template <class F> inline cudaError_t gpu_operand_upload(const gpu_operand<F> &o, F *dev)
{
    if (o.count <= 0)
        return cudaSuccess;
    if (gpu_operand_dense(o))
        return cudaMemcpy(dev, o.host + o.lo, (size_t)o.count * sizeof(F), cudaMemcpyHostToDevice);
    if (o.region == REGION_ALL && o.ld > 0)
        return cudaMemcpy2D(dev - o.lo, (size_t)o.ld * sizeof(F), o.host, (size_t)o.ld * sizeof(F), (size_t)o.rows * sizeof(F), (size_t)o.cols, cudaMemcpyHostToDevice);
    std::vector<F> stage((size_t)o.count);
    gpu_operand_gather(o, stage.data());
    return cudaMemcpy(dev, stage.data(), (size_t)o.count * sizeof(F), cudaMemcpyHostToDevice);
}

template <class F> inline cudaError_t gpu_operand_download(const gpu_operand<F> &o, const F *dev)
{
    if (!o.download || o.count <= 0)
        return cudaSuccess;
    if (gpu_operand_dense(o))
        return cudaMemcpy(const_cast<F *>(o.host) + o.lo, dev, (size_t)o.count * sizeof(F), cudaMemcpyDeviceToHost);
    if (o.region == REGION_ALL && o.ld > 0)
        return cudaMemcpy2D(const_cast<F *>(o.host), (size_t)o.ld * sizeof(F), dev - o.lo, (size_t)o.ld * sizeof(F), (size_t)o.rows * sizeof(F), (size_t)o.cols, cudaMemcpyDeviceToHost);
    std::vector<F> stage((size_t)o.count);
    cudaError_t e = cudaMemcpy(stage.data(), dev, (size_t)o.count * sizeof(F), cudaMemcpyDeviceToHost);
    if (e == cudaSuccess)
        gpu_operand_scatter(o, stage.data());
    return e;
}

} // namespace mplapack_blas
#endif

#endif
