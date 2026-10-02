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
 * Runtime-precision MPFR Rgemm on CUDA (see Rgemm_rt_cuda.h).
 *
 * Ported from mpc_cuda demos/matmul_mpfr.cu: a fixed pool of threads
 * grid-strides over the elements of C; each element is a dot product computed
 * with cu_mpfr_mul/cu_mpfr_add; the limbs that cu_mpfr allocates internally
 * come from a per-thread bump arena that is reset at every step of the dot
 * product.  The accumulators live in a per-thread workspace (the precision
 * is only known at run time, so they cannot be stack arrays as in the demo),
 * and operands are read in place through the MPFR custom interface.
 *
 * The order of operations is that of openmp/Rgemm_*_omp.cpp.
 *
 * Built with -DMPLAPACK_MPFR_CUDA_HOST_EMULATION, gemm_rt runs the same
 * element code on the host instead of launching kernels, with the cu_mpfr
 * calls mapped to the system MPFR (cu_mpfr is MPFR 4.2.2 compiled for the
 * device, and libmpc_cuda's host code needs a GPU for its managed memory).
 * The tests use it to check the element code without a GPU.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Rgemm_rt_cuda.h"

#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include <mpfr.h>
typedef mp_limb_t cu_mp_limb_t;
typedef mpfr_t cu_mpfr_t;
typedef mpfr_ptr cu_mpfr_ptr;
typedef mpfr_srcptr cu_mpfr_srcptr;
typedef mpfr_prec_t cu_mpfr_prec_t;
typedef mpfr_exp_t cu_mpfr_exp_t;
#define CU_MPFR_RNDN MPFR_RNDN
#define CU_MPFR_ZERO_KIND MPFR_ZERO_KIND
#define CU_MPFR_REGULAR_KIND MPFR_REGULAR_KIND
#define cu_mpfr_custom_init_set mpfr_custom_init_set
#define cu_mpfr_custom_get_kind mpfr_custom_get_kind
#define cu_mpfr_custom_get_exp mpfr_custom_get_exp
#define cu_mpfr_custom_get_significand mpfr_custom_get_significand
#define cu_mpfr_mul mpfr_mul
#define cu_mpfr_add mpfr_add
#define cu_mpfr_set mpfr_set
#define cu_mpfr_set_zero mpfr_set_zero
#define RT_HD
#else
#include <cuda_runtime.h>
#include <mpc_cuda.cuh>
#define RT_HD __host__ __device__
#endif

static_assert(sizeof(cu_mp_limb_t) == sizeof(unsigned long long), "64-bit limbs expected");
static_assert(CU_MPFR_ZERO_KIND == mplapack_mpfr_cuda::RT_ZERO_KIND && CU_MPFR_REGULAR_KIND == mplapack_mpfr_cuda::RT_REGULAR_KIND, "MPFR custom kinds differ");

namespace mplapack_mpfr_cuda {

namespace {

// Device-side view of a pack.
struct dpack {
    cu_mp_limb_t *L;
    long *E;
    signed char *K;
};

RT_HD inline void arena_reset()
{
#if defined(__CUDA_ARCH__) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
    mpc_cuda_arena_reset();
#endif
}

// x refers to element idx of P (no copy).
RT_HD inline void view(cu_mpfr_ptr x, long prec, int nl, const dpack &P, size_t idx) { cu_mpfr_custom_init_set(x, (int)P.K[idx], (cu_mpfr_exp_t)P.E[idx], (cu_mpfr_prec_t)prec, P.L + idx * nl); }

// Stores x (zero or regular) into element idx of P.
RT_HD inline void store(cu_mpfr_srcptr x, int nl, dpack &P, size_t idx)
{
    int k = cu_mpfr_custom_get_kind(x);
    P.K[idx] = (signed char)k;
    if (k == CU_MPFR_REGULAR_KIND || k == -CU_MPFR_REGULAR_KIND) {
        P.E[idx] = (long)cu_mpfr_custom_get_exp(x);
        const cu_mp_limb_t *d = (const cu_mp_limb_t *)cu_mpfr_custom_get_significand(x);
        for (int t = 0; t < nl; t++)
            P.L[idx * nl + t] = d[t];
    } else {
        P.E[idx] = 0;
    }
}

RT_HD inline size_t opB_index(const gemm_shape &s, long l, long j) { return s.notb ? (size_t)(l + j * s.ldb) : (size_t)(j + l * s.ldb); }

// AB(l, j) = alpha * op(B)(l, j)   ("temp = alpha * B[...]" in Rgemm_NN/NT_omp)
RT_HD void rt_alpha_opB(const gemm_shape &s, long prec, int nl, const dpack &alpha, const dpack &B, dpack &AB, long l, long j, cu_mp_limb_t *W)
{
    cu_mpfr_t a, b, t;
    cu_mpfr_custom_init_set(t, CU_MPFR_ZERO_KIND, 0, (cu_mpfr_prec_t)prec, W);
    view(a, prec, nl, alpha, 0);
    view(b, prec, nl, B, opB_index(s, l, j));
    arena_reset();
    cu_mpfr_mul(t, a, b, CU_MPFR_RNDN);
    store(t, nl, AB, (size_t)(l + j * s.k));
}

// One element of C; W holds 3*nl limbs of thread-private workspace.
RT_HD void rt_element(const gemm_shape &s, long prec, int nl, const dpack &alpha, const dpack &beta, const dpack &A, const dpack &B, const dpack &AB, dpack &C, long i, long j, cu_mp_limb_t *W)
{
    const cu_mpfr_prec_t p = (cu_mpfr_prec_t)prec;
    cu_mpfr_t c, t, temp, x, y;
    cu_mpfr_custom_init_set(c, CU_MPFR_ZERO_KIND, 0, p, W);
    cu_mpfr_custom_init_set(t, CU_MPFR_ZERO_KIND, 0, p, W + nl);
    cu_mpfr_custom_init_set(temp, CU_MPFR_ZERO_KIND, 0, p, W + 2 * nl);
    const size_t ci = (size_t)(i + j * s.ldc);

    arena_reset();
    if (s.beta_zero) {
        cu_mpfr_set_zero(c, 1); // C[...] = 0.0
    } else {
        view(x, prec, nl, C, ci);
        if (s.beta_one) {
            cu_mpfr_set(c, x, CU_MPFR_RNDN);
        } else {
            view(y, prec, nl, beta, 0);
            cu_mpfr_mul(c, y, x, CU_MPFR_RNDN); // C = beta * C
        }
    }
    if (s.nota) {
        for (long l = 0; l < s.k; l++) {
            arena_reset();
            view(x, prec, nl, AB, (size_t)(l + j * s.k));
            view(y, prec, nl, A, (size_t)(i + l * s.lda));
            cu_mpfr_mul(t, x, y, CU_MPFR_RNDN); // C += temp * A
            cu_mpfr_add(c, c, t, CU_MPFR_RNDN);
        }
    } else {
        cu_mpfr_set_zero(temp, 1); // temp = 0.0
        for (long l = 0; l < s.k; l++) {
            arena_reset();
            view(x, prec, nl, A, (size_t)(l + i * s.lda));
            view(y, prec, nl, B, opB_index(s, l, j));
            cu_mpfr_mul(t, x, y, CU_MPFR_RNDN); // temp += A * B
            cu_mpfr_add(temp, temp, t, CU_MPFR_RNDN);
        }
        arena_reset();
        view(x, prec, nl, alpha, 0);
        cu_mpfr_mul(t, x, temp, CU_MPFR_RNDN); // C += alpha * temp
        cu_mpfr_add(c, c, t, CU_MPFR_RNDN);
    }
    store(c, nl, C, ci);
}

long env_long(const char *name, long fallback)
{
    const char *v = std::getenv(name);
    if (v == NULL || v[0] == '\0')
        return fallback;
    char *end = NULL;
    long r = std::strtol(v, &end, 10);
    return (end != v && *end == '\0' && r > 0) ? r : fallback;
}

#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION

int run(const gemm_shape &s, long prec, int nl, const dpack &alpha, const dpack &beta, const dpack &A, const dpack &B, dpack &AB, dpack &C)
{
    std::vector<cu_mp_limb_t> W(3 * (size_t)nl);
    if (s.nota)
        for (long j = 0; j < s.n; j++)
            for (long l = 0; l < s.k; l++)
                rt_alpha_opB(s, prec, nl, alpha, B, AB, l, j, W.data());
    for (long j = 0; j < s.n; j++)
        for (long i = 0; i < s.m; i++)
            rt_element(s, prec, nl, alpha, beta, A, B, AB, C, i, j, W.data());
    return 0;
}

#else

__global__ void alpha_opB_kernel(gemm_shape s, long prec, int nl, dpack alpha, dpack B, dpack AB, cu_mp_limb_t *work)
{
    const long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    const long stride = (long)gridDim.x * blockDim.x;
    cu_mp_limb_t *W = work + (size_t)tid * 3 * nl;
    for (long e = tid; e < s.k * s.n; e += stride)
        rt_alpha_opB(s, prec, nl, alpha, B, AB, e % s.k, e / s.k, W);
}

__global__ void gemm_kernel(gemm_shape s, long prec, int nl, dpack alpha, dpack beta, dpack A, dpack B, dpack AB, dpack C, cu_mp_limb_t *work)
{
    const long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    const long stride = (long)gridDim.x * blockDim.x;
    cu_mp_limb_t *W = work + (size_t)tid * 3 * nl;
    for (long e = tid; e < s.m * s.n; e += stride)
        rt_element(s, prec, nl, alpha, beta, A, B, AB, C, e % s.m, e / s.m, W);
}

bool verbose()
{
    const char *v = std::getenv("MPLAPACK_MPFR_CUDA_VERBOSE");
    return v != NULL && v[0] != '\0' && v[0] != '0';
}

bool failed(cudaError_t e, const char *what)
{
    if (e == cudaSuccess)
        return false;
    if (verbose())
        std::fprintf(stderr, "mplapack mpfr cuda (runtime precision): %s: %s\n", what, cudaGetErrorString(e));
    return true;
}

// Device copy of a pack of count elements; freed on destruction.
struct device_pack {
    dpack d;
    device_pack() { d.L = NULL, d.E = NULL, d.K = NULL; }
    ~device_pack()
    {
        cudaFree(d.L);
        cudaFree(d.E);
        cudaFree(d.K);
    }
    bool allocate(size_t count, int nl)
    {
        size_t c = count ? count : 1;
        return !failed(cudaMalloc(&d.L, c * nl * sizeof(cu_mp_limb_t)), "cudaMalloc limbs") && !failed(cudaMalloc(&d.E, c * sizeof(long)), "cudaMalloc exponents") && !failed(cudaMalloc(&d.K, c), "cudaMalloc kinds");
    }
    bool upload(const rt_pack &h, size_t count, int nl)
    {
        return !failed(cudaMemcpy(d.L, h.limbs, count * nl * sizeof(cu_mp_limb_t), cudaMemcpyHostToDevice), "copy limbs") && !failed(cudaMemcpy(d.E, h.exp, count * sizeof(long), cudaMemcpyHostToDevice), "copy exponents") && !failed(cudaMemcpy(d.K, h.kind, count, cudaMemcpyHostToDevice), "copy kinds");
    }
    bool download(rt_pack &h, size_t count, int nl) const
    {
        return !failed(cudaMemcpy(h.limbs, d.L, count * nl * sizeof(cu_mp_limb_t), cudaMemcpyDeviceToHost), "copy limbs back") && !failed(cudaMemcpy(h.exp, d.E, count * sizeof(long), cudaMemcpyDeviceToHost), "copy exponents back") && !failed(cudaMemcpy(h.kind, d.K, count, cudaMemcpyDeviceToHost), "copy kinds back");
    }

  private:
    device_pack(const device_pack &);
    device_pack &operator=(const device_pack &);
};

struct device_mem {
    void *p;
    device_mem() : p(NULL) {}
    ~device_mem() { cudaFree(p); }

  private:
    device_mem(const device_mem &);
    device_mem &operator=(const device_mem &);
};

// Raises a device limit for the lifetime of the object and restores it.
struct limit_guard {
    cudaLimit which;
    size_t old;
    bool changed;
    limit_guard(cudaLimit w, size_t want) : which(w), old(0), changed(false)
    {
        if (cudaDeviceGetLimit(&old, w) == cudaSuccess && old < want)
            changed = cudaDeviceSetLimit(w, want) == cudaSuccess;
    }
    ~limit_guard()
    {
        if (changed)
            cudaDeviceSetLimit(which, old);
    }
};

#endif

} // namespace

int gemm_rt(const gemm_shape &s, long prec, const rt_pack &alpha, const rt_pack &beta, const rt_pack &A, size_t countA, const rt_pack &B, size_t countB, rt_pack &C)
{
    const int nl = (int)((prec + 63) / 64);
    const size_t countC = (size_t)s.m * s.n;
    const size_t countAB = s.nota ? (size_t)s.k * s.n : 0;
#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    (void)countA;
    (void)countB;
    std::vector<cu_mp_limb_t> abL(countAB * nl + 1);
    std::vector<long> abE(countAB + 1);
    std::vector<signed char> abK(countAB + 1);
    dpack dal = {(cu_mp_limb_t *)alpha.limbs, alpha.exp, alpha.kind}, dbe = {(cu_mp_limb_t *)beta.limbs, beta.exp, beta.kind};
    dpack dA = {(cu_mp_limb_t *)A.limbs, A.exp, A.kind}, dB = {(cu_mp_limb_t *)B.limbs, B.exp, B.kind};
    dpack dC = {(cu_mp_limb_t *)C.limbs, C.exp, C.kind}, dAB = {abL.data(), abE.data(), abK.data()};
    return run(s, prec, nl, dal, dbe, dA, dB, dAB, dC);
#else
    int ndev = 0;
    if (failed(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount") || ndev == 0)
        return 1;

    // Thread pool and arena as in demos/matmul_mpfr.cu (16 KB per thread at
    // 1024 bits); a thread whose slab overflows falls back to device malloc.
    const long blocks = env_long("MPLAPACK_MPFR_CUDA_RT_BLOCKS", 256);
    const long threads = env_long("MPLAPACK_MPFR_CUDA_RT_THREADS", 32);
    const size_t pool = (size_t)blocks * threads;
    const size_t slab = (size_t)16 * 1024 * ((prec + 1023) / 1024);
    limit_guard stack(cudaLimitStackSize, 64 * 1024);
    limit_guard heap(cudaLimitMallocHeapSize, (size_t)64 * 1024 * 1024);

    device_pack dal, dbe, dA, dB, dC, dAB;
    device_mem arena, top, work;
    if (!dal.allocate(1, nl) || !dbe.allocate(1, nl) || !dA.allocate(countA, nl) || !dB.allocate(countB, nl) || !dC.allocate(countC, nl) || !dAB.allocate(countAB, nl))
        return 2;
    if (failed(cudaMalloc(&arena.p, pool * slab), "cudaMalloc arena") || failed(cudaMalloc(&top.p, pool * sizeof(size_t)), "cudaMalloc arena tops") || failed(cudaMalloc(&work.p, pool * 3 * nl * sizeof(cu_mp_limb_t)), "cudaMalloc workspace"))
        return 2;
    if (failed(cudaMemset(top.p, 0, pool * sizeof(size_t)), "clear arena tops"))
        return 2;
    if (!dal.upload(alpha, 1, nl) || !dbe.upload(beta, 1, nl) || !dA.upload(A, countA, nl) || !dB.upload(B, countB, nl) || (!s.beta_zero && !dC.upload(C, countC, nl)))
        return 3;

    char *saved_base = mpc_cuda_arena_base;
    size_t saved_slab = mpc_cuda_arena_slab;
    size_t *saved_top = mpc_cuda_arena_top;
    mpc_cuda_arena_base = static_cast<char *>(arena.p);
    mpc_cuda_arena_slab = slab;
    mpc_cuda_arena_top = static_cast<size_t *>(top.p);

    cu_mp_limb_t *W = static_cast<cu_mp_limb_t *>(work.p);
    int rc = 0;
    if (s.nota && s.k > 0) {
        alpha_opB_kernel<<<(unsigned)blocks, (unsigned)threads>>>(s, prec, nl, dal.d, dB.d, dAB.d, W);
        if (failed(cudaGetLastError(), "alpha_opB_kernel launch") || failed(cudaDeviceSynchronize(), "alpha_opB_kernel"))
            rc = 4;
    }
    if (rc == 0) {
        gemm_kernel<<<(unsigned)blocks, (unsigned)threads>>>(s, prec, nl, dal.d, dbe.d, dA.d, dB.d, dAB.d, dC.d, W);
        if (failed(cudaGetLastError(), "gemm_kernel launch") || failed(cudaDeviceSynchronize(), "gemm_kernel"))
            rc = 4;
    }
    mpc_cuda_arena_base = saved_base;
    mpc_cuda_arena_slab = saved_slab;
    mpc_cuda_arena_top = saved_top;
    if (rc != 0)
        return rc;
    if (!dC.download(C, countC, nl))
        return 5;
    return 0;
#endif
}

} // namespace mplapack_mpfr_cuda
