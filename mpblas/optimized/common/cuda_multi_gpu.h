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
 * Several CUDA devices for one call whose output columns are independent
 * (GEMM and the like): the columns are cut into one block per device, and
 * one host thread per device computes its block with the single-device code.
 * Every output element is computed exactly as on one device, so the result
 * does not depend on the number of devices.
 *
 * <P>_CUDA_NGPUS=N (N >= 1, <P> = MPLAPACK_DD, MPLAPACK_TD, MPLAPACK_QD,
 * MPLAPACK_MPFR) limits a call to the first N visible devices; by default it
 * uses all of them (CUDA_VISIBLE_DEVICES selects which).  A call with one
 * block runs on the caller's current device, as without this header.
 *
 * Host side only (env_gpus) is plain C++; the rest needs the CUDA runtime.
 */

#ifndef _MPLAPACK_CUDA_MULTI_GPU_H_
#define _MPLAPACK_CUDA_MULTI_GPU_H_

#include <cstdlib>
#include <string>

namespace mplapack_cuda_multi {

// The environment prefix of an element type whose library may use several
// devices (specialized by the libraries); nullptr keeps one device.
template <class Ops> struct prefix_of {
    static constexpr const char *value = nullptr;
};

// The <prefix>_CUDA_NGPUS limit, 0 when unset or invalid (all devices).
inline int env_gpus(const char *prefix) {
    const char *v = std::getenv((std::string(prefix) + "_CUDA_NGPUS").c_str());
    if (v == NULL || v[0] == '\0')
        return 0;
    const long n = std::strtol(v, NULL, 10);
    return n >= 1 && n < 1024 ? (int)n : 0;
}

} // namespace mplapack_cuda_multi

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#include <thread>
#include <vector>

namespace mplapack_cuda_multi {

// Devices a call may use: the visible ones, at most limit (0: no limit).
inline int devices(int limit) {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess || n <= 0)
        return 0;
    return limit > 0 && limit < n ? limit : n;
}

// Number of column blocks for n columns on ndev devices, each at least
// min_cols wide.
inline int blocks(long n, int ndev, long min_cols) {
    long b = min_cols > 0 ? n / min_cols : n;
    if (b > ndev)
        b = ndev;
    return b < 1 ? 1 : (int)b;
}

// Runs fn(j0, j1), returning 0 or an error code, on [0, n) cut into nb
// blocks of a multiple of align columns: in the calling thread when nb == 1,
// otherwise block p in a thread on device p.  Returns the first nonzero code.
template <class Fn> int for_column_blocks(long n, int nb, long align, Fn fn) {
    if (nb <= 1)
        return fn(0L, n);
    long width = (n + nb - 1) / nb;
    if (align > 1)
        width = (width + align - 1) / align * align;
    // The caller may have raised the stack size of its device for deep
    // device code (MPFR); the other devices get at least the same.
    size_t stack = 0;
    if (cudaDeviceGetLimit(&stack, cudaLimitStackSize) != cudaSuccess)
        stack = 0;
    std::vector<int> rc((size_t)nb, 0);
    std::vector<std::thread> workers;
    for (int p = 0; p < nb; p++) {
        const long j0 = p * width, j1 = (p + 1 == nb || j0 + width > n) ? n : j0 + width;
        if (j0 >= j1)
            break;
        workers.emplace_back([&rc, &fn, stack, p, j0, j1]() {
            cudaError_t e = cudaSetDevice(p);
            size_t cur = 0;
            if (e == cudaSuccess && stack > 0 && cudaDeviceGetLimit(&cur, cudaLimitStackSize) == cudaSuccess && cur < stack)
                e = cudaDeviceSetLimit(cudaLimitStackSize, stack);
            rc[(size_t)p] = e != cudaSuccess ? (int)e : fn(j0, j1);
        });
    }
    for (auto &w : workers)
        w.join();
    for (int r : rc)
        if (r != 0)
            return r;
    return 0;
}

// Peer access between all pairs of the first G devices, where supported
// (otherwise the peer copies go through the host).
inline void enable_peer_access(int G) {
    for (int d = 0; d < G; ++d) {
        if (cudaSetDevice(d) != cudaSuccess)
            continue;
        for (int e = 0; e < G; ++e) {
            int can = 0;
            if (e != d && cudaDeviceCanAccessPeer(&can, d, e) == cudaSuccess && can)
                if (cudaDeviceEnablePeerAccess(e, 0) != cudaSuccess)
                    cudaGetLastError(); // already enabled
        }
    }
}

// Restores the calling thread's current device.
struct device_restore {
    int dev = 0;
    device_restore() { cudaGetDevice(&dev); }
    ~device_restore() { cudaSetDevice(dev); }
};

} // namespace mplapack_cuda_multi
#endif

#endif
