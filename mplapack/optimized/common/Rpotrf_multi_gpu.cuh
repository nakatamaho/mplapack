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
 * The blocked Cholesky of Rpotrf_gpu_common.cuh (cholesky()) on several
 * devices.
 *
 * The column blocks of width nb are dealt out cyclically (block b on device
 * b % G).  Lower: every device holds an n x n matrix in which only its own
 * blocks and the panels it receives are used (right() addresses the whole
 * matrix).  Upper: every device keeps its blocks side by side, so that its
 * part of a block row update is one call.  Step b: the owner of block b
 * factors the diagonal block (lower: and solves the panel below it; upper:
 * after the rank update of the diagonal block) and sends the block column to
 * the other devices; then every device updates its blocks right of b (lower:
 * right(); upper: the block row update and solve).  Each device keeps the
 * two-stream lookahead of cholesky(): the next block (on its owner) on the
 * near stream, the others on the far one.
 *
 * Every element goes through the same operations in the same order as in
 * cholesky(), with the same kernels, so the factor and info are those of one
 * device, also when the factorization stops at a non-positive pivot.
 */
#ifndef MPLAPACK_POTRF_MULTI_GPU_CUH
#define MPLAPACK_POTRF_MULTI_GPU_CUH
#include "Rpotrf_gpu_common.cuh"
#include "../../../mpblas/optimized/common/cuda_multi_gpu.h"
#if defined(__CUDACC__)
#include <memory>
#include <thread>
#include <vector>

namespace mplapack_lapack_gpu {

template <class Context> struct chol_device {
    using E = typename Context::E;
    std::unique_ptr<Context> ctx;
    buffer<E> a;
    cudaEvent_t ready{};
    ~chol_device() {
        if (ready)
            cudaEventDestroy(ready);
    }
};

// Upper: the blocks of a device side by side (ld = n), so that its part
// of a block row update is one GEMM and one TRSM; the block columns it
// receives go to two alternating buffers.
template <class Context> struct chol_upper_device {
    using E = typename Context::E;
    std::unique_ptr<Context> ctx;
    buffer<E> local, panel[2];
    cudaEvent_t ready{}, consumed{};
    long cols = 0;
    ~chol_upper_device() {
        if (ready)
            cudaEventDestroy(ready);
        if (consumed)
            cudaEventDestroy(consumed);
    }
};

template <class Context> static int cholesky_upper_multi(int G, long n, long nb, typename Context::E *A, long &info) {
    using E = typename Context::E;
    using namespace mplapack_cuda_multi;
    device_restore restore;
    const long ld = n, blocks = (n + nb - 1) / nb;
    auto owner = [&](long b) { return (int)(b % G); };
    auto width = [&](long b) { return std::min(nb, n - b * nb); };
    auto lcol = [&](long b) { return (b / G) * nb; };
    auto owned_upto = [&](int d, long b) { return b / G + (d <= b % G ? 1 : 0); }; // blocks <= b on d
    enable_peer_access(G);
    std::vector<chol_upper_device<Context>> dev(G);
    for (int d = 0; d < G; ++d) {
        if (int e = cudaSetDevice(d))
            return e;
        auto &x = dev[d];
        for (long b = d; b < blocks; b += G)
            x.cols += width(b);
        x.ctx.reset(new Context(n, nb));
        int e = x.ctx->error;
        if (!e)
            e = x.local.allocate((size_t)ld * std::max(1L, x.cols));
        if (!e)
            e = x.panel[0].allocate((size_t)ld * nb);
        if (!e)
            e = x.panel[1].allocate((size_t)ld * nb);
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
            int e = cudaMemcpy(dev[d].local.data() + lcol(b) * ld, A + b * nb * ld, sizeof(E) * ld * width(b), cudaMemcpyHostToDevice);
            return e ? e : (int)cudaStreamSynchronize(nullptr);
        }))
        return e;
    auto failed = [&]() {
        for (auto &x : dev)
            if (x.ctx->error)
                return x.ctx->error;
        return 0;
    };
    info = 0;
    for (long b = 0; b < blocks; ++b) {
        const long j = b * nb, jb = width(b);
        const int o = owner(b);
        auto &own = dev[o];
        cudaSetDevice(o);
        E *block = own.local.data() + lcol(b) * ld;
        own.ctx->rank(true, jb, j, block, ld, block + j);
        const long i = chol_recursive(*own.ctx, true, jb, block + j, ld);
        if (i) {
            info = j + i;
            break;
        }
        own.ctx->join_update();
        if (!own.ctx->error)
            own.ctx->error = cudaEventRecord(own.consumed, own.ctx->streams[0]);
        // The block column to the other devices, into the buffer each read
        // two steps ago (its last consumed record follows that read).
        for (int d = 0; d < G && !own.ctx->error; ++d) {
            if (d == o)
                continue;
            int e = cudaStreamWaitEvent(own.ctx->streams[0], dev[d].consumed, 0);
            if (!e)
                e = cudaMemcpyPeerAsync(dev[d].panel[b % 2].data(), d, block, o, sizeof(E) * ld * jb, own.ctx->streams[0]);
            if (e)
                own.ctx->error = e;
        }
        if (!own.ctx->error)
            own.ctx->error = cudaEventRecord(own.ready, own.ctx->streams[0]);
        if (int e = failed())
            return e;
        const bool next = b + 1 < blocks;
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
            E *p = d == o ? block : x.panel[b % 2].data(), *local = x.local.data();
            // Blocks > b: the next one first if it is here (near), then the rest.
            const long first = owned_upto(d, b) * nb, cols = x.cols - first;
            const long near = next && owner(b + 1) == d ? width(b + 1) : 0;
            for (int part = 0; part < 2; ++part) {
                const long c = first + (part ? near : 0), w = part ? cols - near : near;
                ctx.far = part != 0;
                if (w > 0) {
                    E *C = local + c * ld;
                    ctx.update(true, jb, w, j, p, C, C + j, ld);
                    ctx.solve(true, jb, w, p + j, ld, C + j);
                }
                if (!part)
                    ctx.release_update();
            }
            ctx.far = false;
        }
        if (int e = failed())
            return e;
    }
    for (int d = 0; d < G; ++d) {
        cudaSetDevice(d);
        dev[d].ctx->finish();
    }
    if (int e = failed())
        return e;
    return each([&](int d, long b) { return (int)cudaMemcpy(A + b * nb * ld, dev[d].local.data() + lcol(b) * ld, sizeof(E) * ld * width(b), cudaMemcpyDeviceToHost); });
}

// A (n x n, ld = n, host staging; the triangle uplo) is overwritten as by
// cholesky(), whose return value info gets.  Returns 0, or a nonzero code
// (A is then undefined).
template <class Context> static int cholesky_multi(int G, bool upper, long n, long nb, typename Context::E *A, long &info) {
    if (upper)
        return cholesky_upper_multi<Context>(G, n, nb, A, info);
    using E = typename Context::E;
    using namespace mplapack_cuda_multi;
    device_restore restore;
    const long ld = n, blocks = (n + nb - 1) / nb;
    auto owner = [&](long b) { return (int)(b % G); };
    auto width = [&](long b) { return std::min(nb, n - b * nb); };
    enable_peer_access(G);
    std::vector<chol_device<Context>> dev(G);
    for (int d = 0; d < G; ++d) {
        if (int e = cudaSetDevice(d))
            return e;
        auto &x = dev[d];
        x.ctx.reset(new Context(n, nb));
        int e = x.ctx->error;
        if (!e)
            e = x.a.allocate((size_t)ld * n);
        if (!e)
            e = cudaEventCreateWithFlags(&x.ready, cudaEventDisableTiming);
        if (e)
            return e;
    }
    // Every device uploads its blocks, at the same time; the other blocks
    // are only ever written by the panels it receives.
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
            int e = cudaMemcpy(dev[d].a.data() + b * nb * ld, A + b * nb * ld, sizeof(E) * ld * width(b), cudaMemcpyHostToDevice);
            // The nonblocking streams must not race the default-stream upload.
            return e ? e : (int)cudaStreamSynchronize(nullptr);
        }))
        return e;
    // The block column of block b from its owner to the other devices.
    auto send = [&](long b) {
        const int o = owner(b);
        auto &own = dev[o];
        for (int d = 0; d < G && !own.ctx->error; ++d)
            if (d != o)
                own.ctx->error = cudaMemcpyPeerAsync(dev[d].a.data() + b * nb * ld, d, own.a.data() + b * nb * ld, o, sizeof(E) * ld * width(b), own.ctx->streams[0]);
        if (!own.ctx->error)
            own.ctx->error = cudaEventRecord(own.ready, own.ctx->streams[0]);
    };
    auto receive = [&](int d, long b) {
        auto &ctx = *dev[d].ctx;
        if (d != owner(b) && !ctx.error)
            ctx.error = cudaStreamWaitEvent(ctx.streams[0], dev[owner(b)].ready, 0);
    };
    auto failed = [&]() {
        for (auto &x : dev)
            if (x.ctx->error)
                return x.ctx->error;
        return 0;
    };
    info = 0;
    for (long b = 0; b < blocks; ++b) {
        const long j = b * nb, jb = width(b), first = j + jb;
        const int o = owner(b);
        auto &own = *dev[o].ctx;
        E *oa = dev[o].a.data();
        cudaSetDevice(o);
        if (upper)
            own.rank(true, jb, j, oa + j * ld, ld, oa + j + j * ld);
        const long i = chol_recursive(own, upper, jb, oa + j + j * ld, ld);
        // The previous steps' updates are complete on every device, as
        // cholesky() finishes before returning a failure.
        if (i) {
            info = j + i;
            break;
        }
        if (!upper && first < n)
            own.solve(false, n - first, jb, oa + j + j * ld, ld, oa + first + j * ld);
        if (upper)
            own.join_update();
        send(b);
        if (int e = failed())
            return e;
        const bool next = b + 1 < blocks;
        for (int d = 0; d < G; ++d) {
            auto &ctx = *dev[d].ctx;
            E *a = dev[d].a.data();
            cudaSetDevice(d);
            if (upper && d != o)
                ctx.join_update();
            receive(d, b);
            // Lower: the far update follows the panel, not the next block's
            // update (cholesky(): join, release, near, far).
            if (!upper) {
                ctx.join_update();
                ctx.release_update();
            }
            // part 0: the next block if it is here (near); part 1: the others.
            for (int part = 0; part < 2; ++part) {
                ctx.far = part != 0;
                for (long c = part ? b + 2 : b + 1; c < blocks; ++c) {
                    if (owner(c) != d || (part && c == b + 1))
                        continue;
                    if (!part && !(next && c == b + 1))
                        break;
                    const long off = c * nb, w = width(c);
                    if (upper) {
                        ctx.update(true, jb, w, j, a + j * ld, a + off * ld, a + j + off * ld, ld);
                        ctx.solve(true, jb, w, a + j + j * ld, ld, a + j + off * ld);
                    } else
                        ctx.right(n, nb, j, jb, off, w, a, ld);
                    if (!part)
                        break;
                }
                if (!part && upper)
                    ctx.release_update();
            }
            ctx.far = false;
        }
        if (int e = failed())
            return e;
    }
    for (int d = 0; d < G; ++d) {
        cudaSetDevice(d);
        dev[d].ctx->finish();
    }
    if (int e = failed())
        return e;
    return each([&](int d, long b) { return (int)cudaMemcpy(A + b * nb * ld, dev[d].a.data() + b * nb * ld, sizeof(E) * ld * width(b), cudaMemcpyDeviceToHost); });
}

} // namespace mplapack_lapack_gpu
#endif
#endif
