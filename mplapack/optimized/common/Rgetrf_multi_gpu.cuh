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
 * The blocked LU of Rgetrf_gpu_common.cuh (lu()) on several devices.
 *
 * The column blocks of width nb are dealt out cyclically: block b lives on
 * device b % G, each device keeping its blocks side by side (ld = m).  Step b:
 * the owner of block b factors it (lu_context::recursive) and sends the
 * panel and its pivots to the other devices; then every device applies the
 * interchanges to its blocks left of b, and the interchanges, the triangular
 * solve and the update to its blocks right of b.  Each device keeps the
 * two-stream lookahead of lu(): the next panel's block (on its owner) is
 * updated on the near stream, the others on the far one.
 *
 * Every element goes through the same operations in the same order as in
 * lu(), with the same kernels and block widths, so the factors, the pivots
 * and info are those of one device.
 */
#ifndef MPLAPACK_GETRF_MULTI_GPU_CUH
#define MPLAPACK_GETRF_MULTI_GPU_CUH
#include "Rgetrf_gpu_common.cuh"
#include "../../../mpblas/optimized/common/cuda_multi_gpu.h"
#if defined(__CUDACC__)
#include <memory>
#include <thread>
#include <vector>

namespace mplapack_lapack_gpu {
using mplapack_cuda_multi::device_restore;
using mplapack_cuda_multi::enable_peer_access;

template <class Blas, class R> struct lu_device {
    using E = typename Blas::E;
    std::unique_ptr<lu_context<Blas, R>> ctx;
    buffer<E> local, panel[2];
    buffer<long> piv;
    cudaEvent_t ready{}, consumed{};
    long cols = 0; // local columns
    ~lu_device() {
        if (ready)
            cudaEventDestroy(ready);
        if (consumed)
            cudaEventDestroy(consumed);
    }
};

// A (m x n, ld = m, host staging) is overwritten by its factors; piv
// (min(m, n)) and info as from lu().  Returns 0, or a nonzero code (A and
// piv are then undefined).
template <class Blas, class R> static int lu_multi(int G, long m, long n, long nb, R sfmin, typename Blas::E *A, long *piv, long &info) {
    using E = typename Blas::E;
    device_restore restore;
    const long mn = std::min(m, n), ld = std::max(1L, m);
    const long blocks = (n + nb - 1) / nb, panels = (mn + nb - 1) / nb;
    auto owner = [&](long b) { return (int)(b % G); };
    auto width = [&](long b) { return std::min(nb, n - b * nb); };
    // Local column of block b on its owner, and the local columns of device d
    // in blocks < b (left) and > b (right).
    auto lcol = [&](long b) { return (b / G) * nb; };
    auto owned_before = [&](int d, long b) { return d < b % G ? (b / G + 1) : b / G; }; // blocks < b on d
    enable_peer_access(G);
    std::vector<lu_device<Blas, R>> dev(G);
    for (int d = 0; d < G; ++d) {
        if (int e = cudaSetDevice(d))
            return e;
        auto &x = dev[d];
        for (long b = d; b < blocks; b += G)
            x.cols += width(b);
        x.ctx.reset(new lu_context<Blas, R>(std::max(m, n), nb, sfmin));
        int e = x.ctx->error;
        if (!e)
            e = x.local.allocate((size_t)ld * std::max(1L, x.cols));
        if (!e)
            e = x.panel[0].allocate((size_t)ld * nb);
        if (!e)
            e = x.panel[1].allocate((size_t)ld * nb);
        if (!e)
            e = x.piv.allocate(std::max(1L, mn));
        if (!e)
            e = cudaEventCreateWithFlags(&x.ready, cudaEventDisableTiming);
        if (!e)
            e = cudaEventCreateWithFlags(&x.consumed, cudaEventDisableTiming);
        if (e)
            return e;
    }
    // Every device uploads its blocks, at the same time.
    std::vector<int> rc(G, 0);
    {
        std::vector<std::thread> up;
        for (int d = 0; d < G; ++d)
            up.emplace_back([&, d]() {
                int e = cudaSetDevice(d);
                for (long b = d; b < blocks && !e; b += G)
                    e = cudaMemcpy(dev[d].local.data() + lcol(b) * ld, A + b * nb * ld, sizeof(E) * ld * width(b), cudaMemcpyHostToDevice);
                // The nonblocking streams must not race the default-stream upload.
                if (!e)
                    e = cudaStreamSynchronize(nullptr);
                rc[d] = e;
            });
        for (auto &t : up)
            t.join();
    }
    for (int d = 0; d < G; ++d) {
        if (rc[d])
            return rc[d];
        cudaSetDevice(d);
        dev[d].ctx->reset_info();
        if (dev[d].ctx->error)
            return dev[d].ctx->error;
    }
    for (long b = 0; b < panels; ++b) {
        const long j = b * nb, jb = std::min(nb, mn - j);
        const int o = owner(b);
        auto &own = dev[o];
        cudaSetDevice(o);
        E *block = own.local.data() + lcol(b) * ld;
        own.ctx->recursive(m - j, jb, block + j, ld, own.piv.data() + j, j);
        // consumed: every read of the panel buffers before this step is done
        // (the far stream is joined; the near stream is in order).
        own.ctx->join_update();
        if (!own.ctx->error)
            own.ctx->error = cudaEventRecord(own.consumed, own.ctx->streams[0]);
        own.ctx->launch(adjust_pivots{own.piv.data() + j, j}, jb);
        // The panel and its pivots to the other devices, on the owner's near
        // stream, into the buffer the receiver last read two steps ago: its
        // last consumed record (step b - 1) follows that read.
        for (int d = 0; d < G && !own.ctx->error; ++d) {
            if (d == o)
                continue;
            auto &x = dev[d];
            int e = cudaStreamWaitEvent(own.ctx->streams[0], x.consumed, 0);
            if (!e)
                e = cudaMemcpyPeerAsync(x.panel[b % 2].data(), d, block, o, sizeof(E) * ld * jb, own.ctx->streams[0]);
            if (!e)
                e = cudaMemcpyPeerAsync(x.piv.data() + j, d, own.piv.data() + j, o, sizeof(long) * jb, own.ctx->streams[0]);
            if (e)
                own.ctx->error = e;
        }
        if (!own.ctx->error)
            own.ctx->error = cudaEventRecord(own.ready, own.ctx->streams[0]);
        for (int d = 0; d < G; ++d) {
            auto &x = dev[d];
            auto &ctx = *x.ctx;
            cudaSetDevice(d);
            if (d != o) {
                ctx.join_update();
                if (!ctx.error)
                    ctx.error = cudaEventRecord(x.consumed, ctx.streams[0]);
                if (!ctx.error)
                    ctx.error = cudaStreamWaitEvent(ctx.streams[0], own.ready, 0);
            }
            if (ctx.error)
                break;
            E *a = d == o ? block : x.panel[b % 2].data(), *local = x.local.data();
            const long *p = x.piv.data();
            // Interchanges left of the panel (blocks < b, all of width nb).
            const long left = owned_before(d, b) * nb;
            if (left)
                ctx.launch(swaps<E>{local, ld, j + 1, j + jb, 1, p}, left);
            // Blocks > b: the next panel's first (near), then the others (far).
            const long first = left + (d == o ? jb : 0), cols = x.cols - first;
            const bool next_here = b + 1 < blocks && owner(b + 1) == d;
            for (int part = 0; part < 2; ++part) {
                const long c = first + (part && next_here ? width(b + 1) : 0);
                const long w = part ? cols - (c - first) : (next_here ? width(b + 1) : 0);
                ctx.far = part != 0;
                if (w > 0) {
                    E *C = local + c * ld;
                    ctx.launch(swaps<E>{C, ld, j + 1, j + jb, 1, p}, w);
                    ctx.solve(jb, w, a + j, C + j, ld);
                    ctx.update(m - j - jb, w, jb, a + j + jb, C + j, C + j + jb, ld);
                }
                if (!part)
                    ctx.release_update();
            }
            ctx.far = false;
        }
        for (auto &x : dev)
            if (x.ctx->error)
                return x.ctx->error;
    }
    // The first singular column: the smallest one recorded on a device.
    info = 0;
    for (int d = 0; d < G; ++d) {
        cudaSetDevice(d);
        const long i = dev[d].ctx->info();
        dev[d].ctx->finish();
        if (dev[d].ctx->error)
            return dev[d].ctx->error;
        if (i && (!info || i < info))
            info = i;
    }
    // All the blocks are final: the pivots (every device has all of them),
    // then every device downloads its blocks, at the same time.
    cudaSetDevice(0);
    if (mn)
        if (int e = cudaMemcpy(piv, dev[0].piv.data(), sizeof(long) * mn, cudaMemcpyDeviceToHost))
            return e;
    {
        std::vector<std::thread> down;
        for (int d = 0; d < G; ++d)
            down.emplace_back([&, d]() {
                int e = cudaSetDevice(d);
                for (long b = d; b < blocks && !e; b += G)
                    e = cudaMemcpy(A + b * nb * ld, dev[d].local.data() + lcol(b) * ld, sizeof(E) * ld * width(b), cudaMemcpyDeviceToHost);
                rc[d] = e;
            });
        for (auto &t : down)
            t.join();
    }
    for (int d = 0; d < G; ++d)
        if (rc[d])
            return rc[d];
    return 0;
}

} // namespace mplapack_lapack_gpu
#endif
#endif
