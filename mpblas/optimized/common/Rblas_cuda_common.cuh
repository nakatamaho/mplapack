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
 * dd/td/qd Level 1/2/3 BLAS on CUDA: the lanes of Rblas_lanes_common.h, one
 * thread per lane, with the arithmetic of Rgemm_cuda_ops_common.h (libQD3's
 * branch-free operations), so that the results are bit-identical to the CPU
 * routines.  Instantiated per type in <type>/cuda/Rgemm_gpu_cuda.cu.  With
 * MPLAPACK_GEMM_CUDA_HOST_EMULATION and a C++ compiler the lanes run on the
 * host instead (tests without a GPU), on copies of the operands made as for
 * the GPU.
 */

#ifndef _MPLAPACK_RBLAS_CUDA_COMMON_CUH_
#define _MPLAPACK_RBLAS_CUDA_COMMON_CUH_

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "Rgemm_cuda_ops_common.h"
#include "Rblas_lanes_common.h"
#include "Rblas_cuda_bridge_common.h"

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#elif !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#error "Rblas_cuda_common.cuh: compile with nvcc, or define MPLAPACK_GEMM_CUDA_HOST_EMULATION"
#endif

namespace mplapack_blas {

// Engine of Rblas_lanes_common.h on the element type of Ops.
template <class Ops> struct ops_engine {
    typedef typename Ops::F V;
    struct M {
        V *p;
        long ld;
    };
    static const bool has_div = true;
    MPLAPACK_BLAS_HD V reg(int) const { return Ops::from_d(0.0); }
    MPLAPACK_BLAS_HD V at(const M &a, long i, long j) const { return a.p[i + j * a.ld]; }
    MPLAPACK_BLAS_HD void put(M &a, long i, long j, const V &v) const { a.p[i + j * a.ld] = v; }
    MPLAPACK_BLAS_HD void set(V &r, const V &a) const { r = a; }
    MPLAPACK_BLAS_HD void zero(V &r) const { r = Ops::from_d(0.0); }
    MPLAPACK_BLAS_HD void one(V &r) const { r = Ops::from_d(1.0); }
    MPLAPACK_BLAS_HD void add(V &r, const V &a, const V &b) const { r = Ops::add(a, b); }
    MPLAPACK_BLAS_HD void sub(V &r, const V &a, const V &b) const { r = Ops::sub(a, b); }
    MPLAPACK_BLAS_HD void mul(V &r, const V &a, const V &b) const { r = Ops::mul(a, b); }
    MPLAPACK_BLAS_HD void div(V &r, const V &a, const V &b) const { r = Ops::div(a, b); }
    // a == 0.0 and a == 1.0 of the C++ types
    MPLAPACK_BLAS_HD bool iszero(const V &a) const { return Ops::eq(a, Ops::from_d(0.0)); }
    MPLAPACK_BLAS_HD bool isone(const V &a) const { return Ops::eq(a, Ops::from_d(1.0)); }
};

template <class Ops> struct gpu_lane_args {
    blas_shape s;
    typename Ops::F *al, *be, *a, *b, *c;
    long lda, ldb, ldc; // of the operands (Rtrmv: b is a copy of c)
};

template <class Ops> MPLAPACK_BLAS_HD void gpu_lane(const gpu_lane_args<Ops> &g, long lane)
{
    typedef ops_engine<Ops> E;
    E e;
    blas_mats<E> x;
    x.al.p = g.al, x.al.ld = 1;
    x.be.p = g.be, x.be.ld = 1;
    x.a.p = g.a, x.a.ld = g.lda;
    x.b.p = g.b, x.b.ld = g.ldb;
    x.c.p = g.c, x.c.ld = g.ldc;
    blas_lane(e, g.s, x, lane);
}

// ---- Rsyrk / Rsyr2k tiled: a TILE x TILE block of threads computes a tile
// of C, the columns l of A (B) in chunks of TILE held in shared memory.  Every
// thread performs the operations of l3_syrk / l3_syr2k (Rblas_lanes_common.h)
// on its element in the same order, so the result is the same bit for bit.
// The phases (load a chunk, use it) are separate functions so that the host
// emulation runs them for all the threads of a block in turn.

static const int BLAS_TILE = 16;

template <class Ops> struct syrk_shared {
    typedef typename Ops::F F;
    // N: x(i, l) and x(j, l) of the chunk; T: x(l, i) and x(l, j)
    F ai[BLAS_TILE][BLAS_TILE + 1], aj[BLAS_TILE][BLAS_TILE + 1]; // [l][t]
    F bi[BLAS_TILE][BLAS_TILE + 1], bj[BLAS_TILE][BLAS_TILE + 1];
};

template <class Ops> struct syrk_state {
    typename Ops::F c, t1, t2;
};

template <class Ops> MPLAPACK_BLAS_HD bool syrk_mine(const blas_shape &s, long i, long j) { return i < s.n && j < s.n && (s.upper ? i <= j : i >= j); }

template <class Ops> MPLAPACK_BLAS_HD void syrk_init(const gpu_lane_args<Ops> &g, syrk_state<Ops> &st, long i, long j)
{
    const blas_shape &s = g.s;
    if (!syrk_mine<Ops>(s, i, j))
        return;
    if (!s.trans) {
        if (s.beta_zero)
            st.c = Ops::from_d(0.0);
        else if (!s.beta_one)
            st.c = Ops::mul(g.be[0], g.c[i + j * g.ldc]);
        else
            st.c = g.c[i + j * g.ldc];
    } else {
        st.t1 = Ops::from_d(0.0);
        st.t2 = Ops::from_d(0.0);
    }
}

// thread (ti, tj) of block (bi, bj) loads element (ti, tj) of each chunk array
template <class Ops> MPLAPACK_BLAS_HD void syrk_load(const gpu_lane_args<Ops> &g, syrk_shared<Ops> &sh, long i0, long j0, long p0, int ti, int tj)
{
    const blas_shape &s = g.s;
    const bool two = s.op == L3_SYR2K;
    const long l = p0 + tj; // the chunk column (N) or row (T) tj
    if (l < s.k) {
        const long ri = i0 + ti, rj = j0 + ti;
        if (ri < s.n) {
            sh.ai[tj][ti] = s.trans ? g.a[l + ri * g.lda] : g.a[ri + l * g.lda];
            if (two)
                sh.bi[tj][ti] = s.trans ? g.b[l + ri * g.ldb] : g.b[ri + l * g.ldb];
        }
        if (rj < s.n) {
            sh.aj[tj][ti] = s.trans ? g.a[l + rj * g.lda] : g.a[rj + l * g.lda];
            if (two)
                sh.bj[tj][ti] = s.trans ? g.b[l + rj * g.ldb] : g.b[rj + l * g.ldb];
        }
    }
}

template <class Ops> MPLAPACK_BLAS_HD void syrk_step(const gpu_lane_args<Ops> &g, const syrk_shared<Ops> &sh, syrk_state<Ops> &st, long i, long j, long kc, int ti, int tj)
{
    typedef typename Ops::F F;
    const blas_shape &s = g.s;
    if (!syrk_mine<Ops>(s, i, j))
        return;
    const F al = g.al[0];
    const F zero = Ops::from_d(0.0);
    for (long l = 0; l < kc; l++) {
        const F &ail = sh.ai[l][ti], &ajl = sh.aj[l][tj];
        if (s.op == L3_SYRK) {
            if (!s.trans) {
                if (!Ops::eq(ajl, zero)) {
                    const F t = Ops::mul(al, ajl);         // temp = alpha * a(j,l)
                    st.c = Ops::add(st.c, Ops::mul(t, ail)); // c(i,j) += temp * a(i,l)
                }
            } else {
                st.t1 = Ops::add(st.t1, Ops::mul(ail, ajl)); // temp += a(l,i) * a(l,j)
            }
        } else {
            const F &bil = sh.bi[l][ti], &bjl = sh.bj[l][tj];
            if (!s.trans) {
                if (!Ops::eq(ajl, zero) || !Ops::eq(bjl, zero)) {
                    const F t1 = Ops::mul(al, bjl), t2 = Ops::mul(al, ajl);
                    st.c = Ops::add(st.c, Ops::add(Ops::mul(ail, t1), Ops::mul(bil, t2)));
                }
            } else {
                st.t1 = Ops::add(st.t1, Ops::mul(ail, bjl)); // temp1 += a(l,i) * b(l,j)
                st.t2 = Ops::add(st.t2, Ops::mul(bil, ajl)); // temp2 += b(l,i) * a(l,j)
            }
        }
    }
}

template <class Ops> MPLAPACK_BLAS_HD void syrk_finish(const gpu_lane_args<Ops> &g, syrk_state<Ops> &st, long i, long j)
{
    typedef typename Ops::F F;
    const blas_shape &s = g.s;
    if (!syrk_mine<Ops>(s, i, j))
        return;
    F c = st.c;
    if (s.trans) {
        const F al = g.al[0], be = g.be[0];
        if (s.op == L3_SYRK) {
            c = Ops::mul(al, st.t1); // c = alpha*temp [+ beta*c]
            if (!s.beta_zero)
                c = Ops::add(c, Ops::mul(be, g.c[i + j * g.ldc]));
        } else if (s.beta_zero) {
            c = Ops::add(Ops::mul(al, st.t1), Ops::mul(al, st.t2));
        } else {
            c = Ops::mul(be, g.c[i + j * g.ldc]);
            c = Ops::add(c, Ops::mul(al, st.t1));
            c = Ops::add(c, Ops::mul(al, st.t2));
        }
    }
    g.c[i + j * g.ldc] = c;
}

// the tiled form is used for Rsyrk and Rsyr2k
template <class Ops> inline bool syrk_tiled(const blas_shape &s) { return s.op == L3_SYRK || s.op == L3_SYR2K; }

#if defined(__CUDACC__)

template <class Ops> __global__ void k_blas_lanes(gpu_lane_args<Ops> g, long lanes)
{
    const long lane = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (lane < lanes)
        gpu_lane<Ops>(g, lane);
}

template <class Ops> __global__ void k_syrk_tiled(gpu_lane_args<Ops> g)
{
    __shared__ syrk_shared<Ops> sh;
    const int ti = threadIdx.x, tj = threadIdx.y;
    const long i0 = (long)blockIdx.x * BLAS_TILE, j0 = (long)blockIdx.y * BLAS_TILE;
    // blocks entirely outside the referenced triangle have nothing to do
    if (g.s.upper ? i0 > j0 + BLAS_TILE - 1 : j0 > i0 + BLAS_TILE - 1)
        return;
    const long i = i0 + ti, j = j0 + tj;
    syrk_state<Ops> st;
    syrk_init<Ops>(g, st, i, j);
    for (long p0 = 0; p0 < g.s.k; p0 += BLAS_TILE) {
        syrk_load<Ops>(g, sh, i0, j0, p0, ti, tj);
        __syncthreads();
        const long kc = g.s.k - p0 < BLAS_TILE ? g.s.k - p0 : BLAS_TILE;
        syrk_step<Ops>(g, sh, st, i, j, kc, ti, tj);
        __syncthreads();
    }
    syrk_finish<Ops>(g, st, i, j);
}

// Device copy of the span of an operand (only its referenced elements are
// transferred, Rblas_gpu_operand_common.h).
template <class F> class gpu_span {
  public:
    gpu_span() : d_(NULL) {}
    ~gpu_span() { cudaFree(d_); }
    cudaError_t up(const gpu_operand<F> &o)
    {
        o_ = o;
        cudaError_t e = cudaMalloc(&d_, (o.count ? o.count : 1) * sizeof(F));
        if (e == cudaSuccess)
            e = gpu_operand_upload(o, d_);
        return e;
    }
    cudaError_t down() { return gpu_operand_download(o_, d_); }
    F *at0() const { return d_ - o_.lo; } // logical element (0, 0)

  private:
    F *d_;
    gpu_operand<F> o_;
    gpu_span(const gpu_span &);
    gpu_span &operator=(const gpu_span &);
};

template <class Ops> int blas_gpu(const blas_shape &s, const typename Ops::F &alpha, const typename Ops::F &beta, const gpu_operand<typename Ops::F> &a, const gpu_operand<typename Ops::F> &b, const gpu_operand<typename Ops::F> &c)
{
    typedef typename Ops::F F;
    gpu_operand<F> sc[2];
    sc[0] = gpu_scalar_operand(&alpha);
    sc[1] = gpu_scalar_operand(&beta);
    gpu_span<F> dal, dbe, da, db, dc;
    cudaError_t e = dal.up(sc[0]);
    if (e == cudaSuccess)
        e = dbe.up(sc[1]);
    if (e == cudaSuccess)
        e = da.up(a);
    if (e == cudaSuccess)
        e = db.up(b);
    if (e == cudaSuccess)
        e = dc.up(c);
    if (e != cudaSuccess)
        return (int)e;
    gpu_lane_args<Ops> g;
    g.s = s;
    g.al = dal.at0(), g.be = dbe.at0(), g.a = da.at0(), g.b = db.at0(), g.c = dc.at0();
    g.lda = a.ld, g.ldb = b.ld, g.ldc = c.ld;
    const long lanes = blas_lanes(s);
    const int threads = 64;
    if (syrk_tiled<Ops>(s)) {
        const unsigned nb = (unsigned)((s.n + BLAS_TILE - 1) / BLAS_TILE);
        k_syrk_tiled<Ops><<<dim3(nb, nb), dim3(BLAS_TILE, BLAS_TILE)>>>(g);
        e = cudaGetLastError();
        if (e == cudaSuccess)
            e = cudaDeviceSynchronize();
    } else if (lanes > 0) {
        k_blas_lanes<Ops><<<(unsigned)((lanes + threads - 1) / threads), threads>>>(g, lanes);
        e = cudaGetLastError();
        if (e == cudaSuccess)
            e = cudaDeviceSynchronize();
    }
    if (e == cudaSuccess)
        e = dc.down();
    if (e == cudaSuccess)
        e = da.down();
    return (int)e;
}

#else

// Host emulation: the same lanes on copies of the operand spans.
template <class F> struct emul_span {
    std::vector<F> v;
    gpu_operand<F> o;
    void up(const gpu_operand<F> &op)
    {
        o = op;
        v.assign((size_t)(op.count ? op.count : 1), F());
        if (op.count)
            gpu_operand_gather(op, v.data());
    }
    void down()
    {
        if (o.download && o.count)
            gpu_operand_scatter(o, v.data());
    }
    F *at0() { return v.data() - o.lo; }
};

template <class Ops> int blas_gpu(const blas_shape &s, const typename Ops::F &alpha, const typename Ops::F &beta, const gpu_operand<typename Ops::F> &a, const gpu_operand<typename Ops::F> &b, const gpu_operand<typename Ops::F> &c)
{
    typedef typename Ops::F F;
    emul_span<F> dal, dbe, da, db, dc;
    dal.up(gpu_scalar_operand(&alpha));
    dbe.up(gpu_scalar_operand(&beta));
    da.up(a), db.up(b), dc.up(c);
    gpu_lane_args<Ops> g;
    g.s = s;
    g.al = dal.at0(), g.be = dbe.at0(), g.a = da.at0(), g.b = db.at0(), g.c = dc.at0();
    g.lda = a.ld, g.ldb = b.ld, g.ldc = c.ld;
    const long lanes = blas_lanes(s);
    if (syrk_tiled<Ops>(s)) {
        // the tiled kernel, the threads of a block phase by phase
        const long nb = (s.n + BLAS_TILE - 1) / BLAS_TILE;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (long b = 0; b < nb * nb; b++) {
            const long i0 = (b % nb) * BLAS_TILE, j0 = (b / nb) * BLAS_TILE;
            if (s.upper ? i0 > j0 + BLAS_TILE - 1 : j0 > i0 + BLAS_TILE - 1)
                continue;
            std::vector<syrk_shared<Ops> > shv(1);
            syrk_shared<Ops> &sh = shv[0];
            std::vector<syrk_state<Ops> > st((size_t)(BLAS_TILE * BLAS_TILE));
            for (int tj = 0; tj < BLAS_TILE; tj++)
                for (int ti = 0; ti < BLAS_TILE; ti++)
                    syrk_init<Ops>(g, st[(size_t)(ti + tj * BLAS_TILE)], i0 + ti, j0 + tj);
            for (long p0 = 0; p0 < s.k; p0 += BLAS_TILE) {
                for (int tj = 0; tj < BLAS_TILE; tj++)
                    for (int ti = 0; ti < BLAS_TILE; ti++)
                        syrk_load<Ops>(g, sh, i0, j0, p0, ti, tj);
                const long kc = s.k - p0 < BLAS_TILE ? s.k - p0 : BLAS_TILE;
                for (int tj = 0; tj < BLAS_TILE; tj++)
                    for (int ti = 0; ti < BLAS_TILE; ti++)
                        syrk_step<Ops>(g, sh, st[(size_t)(ti + tj * BLAS_TILE)], i0 + ti, j0 + tj, kc, ti, tj);
            }
            for (int tj = 0; tj < BLAS_TILE; tj++)
                for (int ti = 0; ti < BLAS_TILE; ti++)
                    syrk_finish<Ops>(g, st[(size_t)(ti + tj * BLAS_TILE)], i0 + ti, j0 + tj);
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (long lane = 0; lane < lanes; lane++)
            gpu_lane<Ops>(g, lane);
    }
    dc.down();
    da.down();
    return 0;
}

#endif

} // namespace mplapack_blas

#endif
