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
 * factor_context::lu() and factor_context::chol() (Rfactor_device_body.cuh)
 * on several devices, for the fixed-precision real engine (the runtime and
 * complex engines index process-wide __managed__ arenas and stay on one
 * device).  The scheme of the dd/td/qd ones
 * (mplapack/optimized/common/Rgetrf_multi_gpu.cuh, Rpotrf_multi_gpu.cuh):
 *
 * The column blocks of width nb are dealt out cyclically (block b on device
 * b % G).  LU and upper Cholesky: every device keeps its blocks side by side
 * (ld = m); the owner of block b factors it (LU: the panel; upper: the rank
 * update and the diagonal block) and sends the block column to the other
 * devices, into one of two alternating buffers; every device then updates
 * its blocks right of b, in one call per stream (LU: also the interchanges
 * of its blocks left of b).  Lower Cholesky: right() addresses the whole
 * matrix, so every device holds all of it, using only its blocks and the
 * block columns it receives (the owner factors the diagonal block and
 * solves the panel below it), and updates its blocks one by one.  Each
 * device keeps the two-stream lookahead of factor_context.
 *
 * Every element goes through the same operations in the same order as on
 * one device, with the same kernels, so the factors, the pivots and info
 * are those of factor_context::factor().
 */
#ifndef MPLAPACK_MPFR_FACTOR_MULTI_GPU_CUH
#define MPLAPACK_MPFR_FACTOR_MULTI_GPU_CUH
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include "Rfactor_device_body.cuh"
#include "../../common/cuda_multi_gpu.h"
#include <memory>
#include <thread>
#include <vector>

namespace mplapack_mpfr_cuda {
namespace {

template <class E> struct factor_device {
    using V = typename E::V;
    std::unique_ptr<factor_stack_guard> stack;
    factor_buffer<V> constants, local, panel[2];
    std::unique_ptr<factor_context<E>> ctx;
    cudaEvent_t ready{}, consumed{};
    long cols = 0;
    ~factor_device() {
        if (ready)
            cudaEventDestroy(ready);
        if (consumed)
            cudaEventDestroy(consumed);
    }
};

// a: the m x n matrix (ld = m, host), overwritten as by factor(); piv
// (LU, min(m, n)) and info likewise.  stack_bytes: the device stack the
// engine needs (raised on every device).  Returns 0, or a nonzero code (a and piv
// are then undefined).
template <class E> static int factor_multi(int G, const factor_shape &s, typename E::V *a, const typename E::V &safe, long *piv, long &info, size_t stack_bytes) {
    using V = typename E::V;
    using M = typename E::M;
    using namespace mplapack_cuda_multi;
    device_restore restore;
    const long m = s.m, n = s.n, nb = s.nb, mn = std::min(m, n), ld = m;
    const long blocks = (n + nb - 1) / nb;
    const bool whole = s.chol && !s.upper; // lower Cholesky: the whole matrix
    auto owner = [&](long b) { return (int)(b % G); };
    auto width = [&](long b) { return std::min(nb, n - b * nb); };
    // Column of block b in its device's storage.
    auto lcol = [&](long b) { return whole ? b * nb : (b / G) * nb; };
    auto owned_before = [&](int d, long b) { return b / G + (d < b % G ? 1 : 0); }; // blocks < b on d
    enable_peer_access(G);
    std::vector<factor_device<E>> dev(G);
    for (int d = 0; d < G; ++d) {
        if (int e = cudaSetDevice(d))
            return e;
        auto &x = dev[d];
        for (long b = d; b < blocks; b += G)
            x.cols += width(b);
        x.stack.reset(new factor_stack_guard(stack_bytes));
        int e = x.stack->error;
        V c[4] = {V::from_double(1), V::from_double(-1), safe, V::from_double(0)};
        if (!e)
            e = x.constants.allocate(4);
        if (!e)
            e = x.constants.upload(c, 4);
        if (!e)
            e = x.local.allocate((size_t)ld * (whole ? n : std::max(1L, x.cols)));
        if (!e && !whole)
            e = x.panel[0].allocate((size_t)ld * nb);
        if (!e && !whole)
            e = x.panel[1].allocate((size_t)ld * nb);
        if (!e)
            x.ctx.reset(new factor_context<E>(E{}, M{x.constants.p, 4}, m, n));
        if (!e)
            e = x.ctx->error;
        if (!e)
            e = cudaEventCreateWithFlags(&x.ready, cudaEventDisableTiming);
        if (!e)
            e = cudaEventCreateWithFlags(&x.consumed, cudaEventDisableTiming);
        if (e)
            return e;
    }
    std::vector<int> rc(G, 0);
    auto each = [&](auto f) {
        std::vector<std::thread> t;
        for (int d = 0; d < G; ++d)
            t.emplace_back([&, d]() {
                int e = cudaSetDevice(d);
                for (long b = d; b < blocks && !e; b += G)
                    e = f(d, b);
                rc[d] = e;
            });
        for (auto &x : t)
            x.join();
        for (int d = 0; d < G; ++d)
            if (rc[d])
                return rc[d];
        return 0;
    };
    if (int e = each([&](int d, long b) {
            int e = cudaMemcpy(dev[d].local.p + lcol(b) * ld, a + b * nb * ld, sizeof(V) * ld * width(b), cudaMemcpyHostToDevice);
            // The nonblocking streams must not race the default-stream upload.
            return e ? e : (int)cudaStreamSynchronize(nullptr);
        }))
        return e;
    auto failed = [&]() {
        for (auto &x : dev)
            if (x.ctx->error)
                return x.ctx->error;
        return 0;
    };
    // The block column of block b, from its owner's storage to the other
    // devices (whole: at its place; else into the buffer each last read two
    // steps ago, after its consumed record).
    auto send = [&](long b) {
        const int o = owner(b);
        auto &own = *dev[o].ctx;
        V *src = dev[o].local.p + lcol(b) * ld;
        for (int d = 0; d < G && !own.error; ++d) {
            if (d == o)
                continue;
            int e = whole ? 0 : cudaStreamWaitEvent(own.streams[0], dev[d].consumed, 0);
            if (!e)
                e = cudaMemcpyPeerAsync(whole ? dev[d].local.p + lcol(b) * ld : dev[d].panel[b % 2].p, d, src, o, sizeof(V) * ld * width(b), own.streams[0]);
            if (!e && !s.chol && width(b))
                e = cudaMemcpyPeerAsync(dev[d].ctx->pivots.p + b * nb, d, own.pivots.p + b * nb, o, sizeof(long) * std::min(nb, mn - b * nb), own.streams[0]);
            if (e)
                own.error = e;
        }
        if (!own.error)
            own.error = cudaEventRecord(dev[o].ready, own.streams[0]);
    };
    // The receiver's near stream: its earlier reads of the buffers are done
    // (consumed) and the block column has arrived.
    auto receive = [&](int d, long b) {
        auto &ctx = *dev[d].ctx;
        if (d == owner(b) || ctx.error)
            return;
        ctx.join();
        if (!ctx.error)
            ctx.error = cudaEventRecord(dev[d].consumed, ctx.streams[0]);
        if (!ctx.error)
            ctx.error = cudaStreamWaitEvent(ctx.streams[0], dev[owner(b)].ready, 0);
    };
    info = 0;
    const long steps = s.chol ? blocks : (mn + nb - 1) / nb;
    for (long b = 0; b < steps; ++b) {
        const long j = b * nb, jb = s.chol ? width(b) : std::min(nb, mn - j);
        const int o = owner(b);
        auto &own = *dev[o].ctx;
        cudaSetDevice(o);
        M block{dev[o].local.p + lcol(b) * ld, ld};
        if (!s.chol) {
            const long i = own.lu_recursive(m - j, jb, E::shift(block, j), own.pivots.p + j);
            if (!info && i)
                info = j + i;
            own.join();
            if (!own.error)
                own.error = cudaEventRecord(dev[o].consumed, own.streams[0]);
            own.adjust(own.pivots.p + j, jb, j);
        } else {
            if (s.upper)
                own.rank(true, jb, j, block, E::shift(block, j));
            const long i = own.chol_recursive(s.upper, jb, E::shift(block, j));
            // As factor(): the previous steps' updates complete, then stop.
            if (own.error || i) {
                info = i ? j + i : 0;
                break;
            }
            if (!s.upper && j + jb < n)
                own.solve(false, n - j - jb, jb, E::shift(block, j), E::shift(block, j + jb), true);
            if (s.upper) {
                own.join();
                if (!own.error)
                    own.error = cudaEventRecord(dev[o].consumed, own.streams[0]);
            }
        }
        send(b);
        if (int e = failed())
            return e;
        const bool next_here_any = b + 1 < blocks;
        for (int d = 0; d < G; ++d) {
            auto &x = dev[d];
            auto &ctx = *x.ctx;
            cudaSetDevice(d);
            receive(d, b);
            M p = d == o || whole ? M{x.local.p + lcol(b) * ld, ld} : M{x.panel[b % 2].p, ld};
            const bool next_here = next_here_any && owner(b + 1) == d;
            if (whole) {
                // Lower Cholesky (factor_context::chol): join, release, the
                // next block (near), the others (far), one by one.
                if (d == o)
                    ctx.join();
                ctx.release();
                for (long c = b + 1; c < blocks; ++c) {
                    if (owner(c) != d)
                        continue;
                    ctx.far = c != b + 1;
                    ctx.right(n, nb, j, jb, c * nb, width(c), M{x.local.p, ld});
                }
                ctx.far = false;
                continue;
            }
            const long left = owned_before(d, b) * nb;
            const long first = left + (d == o ? (s.chol ? width(b) : jb) : 0), cols = x.cols - first;
            if (!s.chol && left)
                ctx.swap(left, M{x.local.p, ld}, j, j + jb, ctx.pivots.p);
            const long near = next_here ? width(b + 1) : 0;
            for (int part = 0; part < 2; ++part) {
                const long c = first + (part ? near : 0), w = part ? cols - near : near;
                ctx.far = part != 0;
                if (w > 0) {
                    M C{x.local.p + c * ld, ld};
                    if (!s.chol) {
                        ctx.swap(w, C, j, j + jb, ctx.pivots.p);
                        M right = E::shift(C, j);
                        ctx.solve(false, jb, w, E::shift(p, j), right, false);
                        auto job = ctx.job(1, E::shift(p, j + jb));
                        job.b = right;
                        job.c = E::shift(right, jb);
                        job.g = {1, 1, m - j - jb, w, jb, ld, ld, ld, 0, 1};
                        ctx.launch(job, (m - j - jb) * w);
                    } else {
                        ctx.update(true, jb, w, j, p, C, E::shift(C, j));
                        ctx.solve(true, jb, w, E::shift(p, j), E::shift(C, j), true);
                    }
                }
                if (!part)
                    ctx.release();
            }
            ctx.far = false;
        }
        if (int e = failed())
            return e;
    }
    int flag = 0;
    for (int d = 0; d < G; ++d) {
        cudaSetDevice(d);
        auto &ctx = *dev[d].ctx;
        ctx.finish();
        int f = 0;
        if (!ctx.error)
            ctx.error = ctx.failed.download(&f, 1);
        flag |= f;
    }
    if (int e = failed())
        return e;
    if (flag)
        return 1;
    if (!s.chol && mn) {
        cudaSetDevice(0);
        if (dev[0].ctx->pivots.download(piv, mn))
            return 1;
    }
    int e = each([&](int d, long b) { return (int)cudaMemcpy(a + b * nb * ld, dev[d].local.p + lcol(b) * ld, sizeof(V) * ld * width(b), cudaMemcpyDeviceToHost); });
    for (int d = 0; d < G; ++d) {
        cudaSetDevice(d);
        dev[d].stack->restore();
    }
    return e;
}

} // namespace
} // namespace mplapack_mpfr_cuda
#endif
#endif
