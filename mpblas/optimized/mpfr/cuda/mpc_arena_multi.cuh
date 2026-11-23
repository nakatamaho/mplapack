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
 * libmpc_cuda's arena (mpc_cuda_arena_base, _slab, _top: process-wide
 * __managed__ variables; the slot of a thread is its physical index
 * (bx + gx*by + gx*gy*bz)*bsz + tib) shared by several devices: one managed
 * arena of G * P slabs, the P slabs and offsets of device d resident on d.
 * The kernels of device d start d * P threads late (their first d * P
 * threads return at once), so the devices use disjoint slabs.  Needs
 * concurrent managed access on every device; allocate() returns -1
 * otherwise (nothing allocated: run on one device).
 */

#ifndef MPLAPACK_MPFR_MPC_ARENA_MULTI_CUH
#define MPLAPACK_MPFR_MPC_ARENA_MULTI_CUH

#include <cuda_runtime.h>
#include <cstddef>

namespace mplapack_mpfr_cuda {

struct mpc_arena_multi {
    char *base = nullptr;
    size_t *top = nullptr;
    mpc_arena_multi() = default;
    mpc_arena_multi(const mpc_arena_multi &) = delete;
    mpc_arena_multi &operator=(const mpc_arena_multi &) = delete;
    ~mpc_arena_multi() {
        cudaFree(base);
        cudaFree(top);
    }
    // Leaves the current device as it was.
    int allocate(int G, size_t P, size_t slab) {
        int cur = 0;
        if (int e = cudaGetDevice(&cur))
            return e;
        const int e = place(G, P, slab);
        const int r = cudaSetDevice(cur);
        return e ? e : r;
    }

  private:
    int place(int G, size_t P, size_t slab) {
        for (int d = 0; d < G; ++d) {
            int managed = 0;
            if (cudaDeviceGetAttribute(&managed, cudaDevAttrConcurrentManagedAccess, d) || !managed)
                return -1;
        }
        if (int e = cudaMallocManaged(&base, G * P * slab))
            return e;
        if (int e = cudaMallocManaged(&top, G * P * sizeof(size_t)))
            return e;
        for (int d = 0; d < G; ++d) {
            char *b = base + d * P * slab;
            size_t *t = top + d * P;
            int e = cudaSetDevice(d);
            if (!e)
                e = cudaMemAdvise(b, P * slab, cudaMemAdviseSetPreferredLocation, d);
            if (!e)
                e = cudaMemAdvise(t, P * sizeof(size_t), cudaMemAdviseSetPreferredLocation, d);
            if (!e)
                e = cudaMemPrefetchAsync(b, P * slab, d, nullptr);
            if (!e)
                e = cudaMemsetAsync(t, 0, P * sizeof(size_t), nullptr);
            if (!e)
                e = cudaDeviceSynchronize();
            if (e)
                return e;
        }
        return 0;
    }
};

} // namespace mplapack_mpfr_cuda

#endif
