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
 * Test of the dd/qd CUDA Rgemm (compile with -DRGEMM_TEST_DD or
 * -DRGEMM_TEST_QD, together with <type>/cuda/Rgemm_bridge_cuda.cpp and
 * <type>/cuda/Rgemm_gpu_cuda.cu, and link with libmplapack_<type>_opt).
 * Rgemm_gpu() is compared with the CPU Rgemm of libmplapack_<type>_opt:
 *
 *  - naive and tiled kernels: bit-identical to the conventional CPU Rgemm
 *    (k == 0 is never sent to the GPU: m*n*k is below any threshold);
 *  - Ozaki: bit-identical to the CPU Ozaki path (<P>_GEMM_OZAKI=1);
 *  - Winograd (cutoffs 1, 2, 3, 8): exact on small integer matrices, within
 *    2^(20-p) (k |alpha| max|A| max|B| + |beta| max|C|) on random reals;
 *
 * for all transpose combinations, odd sizes, leading dimensions larger than
 * needed and alpha/beta in {general, 0, 1}.  With
 * MPLAPACK_GEMM_CUDA_HOST_EMULATION the GPU code runs on the host; otherwise
 * the test is skipped (exit 77) when no device is present.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <vector>

#if defined(RGEMM_TEST_DD)
#include <mpblas_dd.h>
typedef dd_real T;
static const char *PREFIX = "MPLAPACK_DD";
static const int PREC = 106;
#elif defined(RGEMM_TEST_QD)
#include <mpblas_qd.h>
typedef qd_real T;
static const char *PREFIX = "MPLAPACK_QD";
static const int PREC = 212;
#else
#error "define RGEMM_TEST_DD or RGEMM_TEST_QD"
#endif

bool Rgemm_gpu(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, T alpha, T *A, mplapackint lda, T *B, mplapackint ldb, T beta, T *C, mplapackint ldc);

#ifndef MPLAPACK_GEMM_CUDA_HOST_EMULATION
#include <cuda_runtime.h>
#endif

static std::mt19937_64 rng(20261005);
static double u11() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }
static double todbl(const T &x) { return x.x[0]; }

static T random_real() {
    T x = (T)u11();
    x += (T)(u11() * std::ldexp(1.0, -53));
    x += (T)(u11() * std::ldexp(1.0, -106));
    x += (T)(u11() * std::ldexp(1.0, -159));
    if (rng() % 17 == 0)
        x = (T)0.0;
    return x;
}
static T random_int() { return (T)(double)((long)(rng() % 129) - 64); }

static void env(const char *name, const char *value) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%s%s", PREFIX, name);
    if (value)
        setenv(buf, value, 1);
    else
        unsetenv(buf);
}

static void reset_env() {
    const char *names[] = {"_CUDA", "_CUDA_KERNEL", "_CUDA_WINOGRAD_CUTOFF", "_CUDA_OZAKI", "_CUDA_OZAKI_SLICES", "_GEMM_OZAKI", "_GEMM_WINOGRAD_CUTOFF", "_GEMM_BLOCKED"};
    for (const char *nm : names)
        env(nm, NULL);
    env("_CUDA_MIN_MNK", "1");
    env("_GEMM_MIN_MNK", "1");
    env("_CUDA_VERBOSE", "1");
}

struct problem {
    char ta, tb;
    long m, n, k, lda, ldb, ldc;
    T alpha, beta;
    std::vector<T> A, B, C;
};

static problem make(char ta, char tb, long m, long n, long k, int betamode, bool integers) {
    problem p;
    p.ta = ta;
    p.tb = tb;
    p.m = m;
    p.n = n;
    p.k = k;
    const long ra = ta == 'n' ? m : k, ca = ta == 'n' ? k : m;
    const long rb = tb == 'n' ? k : n, cb = tb == 'n' ? n : k;
    p.lda = ra + 3;
    p.ldb = rb + 1;
    p.ldc = m + 2;
    p.A.resize(p.lda * ca + 1);
    p.B.resize(p.ldb * cb + 1);
    p.C.resize(p.ldc * n + 1);
    for (auto &x : p.A)
        x = integers ? random_int() : random_real();
    for (auto &x : p.B)
        x = integers ? random_int() : random_real();
    for (auto &x : p.C)
        x = integers ? random_int() : random_real();
    p.alpha = integers ? (T)3.0 : random_real();
    if (p.alpha == 0.0)
        p.alpha = (T)0.75;
    p.beta = betamode == 0 ? (T)0.0 : betamode == 1 ? (T)1.0 : (integers ? (T)-2.0 : random_real());
    return p;
}

static std::vector<T> cpu(const problem &p) {
    problem q = p;
    const char ta[2] = {p.ta, 0}, tb[2] = {p.tb, 0};
    Rgemm(ta, tb, q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc);
    return q.C;
}

static long failures = 0, checks = 0;

static std::vector<T> gpu(const problem &p, const char *what) {
    problem q = p;
    checks++;
    if (!Rgemm_gpu(p.ta == 'n', p.tb == 'n', q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc)) {
        failures++;
        std::printf("FAIL %s: %c%c m=%ld n=%ld k=%ld did not run on the GPU\n", what, p.ta, p.tb, p.m, p.n, p.k);
    }
    return q.C;
}

static void expect_identical(const char *what, const problem &p, const std::vector<T> &got, const std::vector<T> &ref) {
    checks++;
    if (std::memcmp(got.data(), ref.data(), sizeof(T) * got.size()) != 0) {
        failures++;
        std::printf("FAIL %s: %c%c m=%ld n=%ld k=%ld not bit-identical\n", what, p.ta, p.tb, p.m, p.n, p.k);
    }
}

static double maxabs(const std::vector<T> &v) {
    double r = 0.0;
    for (const auto &x : v)
        r = std::fmax(r, std::fabs(todbl(x)));
    return r;
}

static void expect_close(const char *what, const problem &p, const std::vector<T> &got, const std::vector<T> &ref) {
    checks++;
    const double bound = std::ldexp(1.0, 20 - PREC) * (p.k * std::fabs(todbl(p.alpha)) * maxabs(p.A) * maxabs(p.B) + std::fabs(todbl(p.beta)) * maxabs(p.C));
    double worst = 0.0;
    for (size_t e = 0; e < got.size(); e++) {
        const long i = (long)(e % p.ldc), j = (long)(e / p.ldc);
        if (i < p.m && j < p.n)
            worst = std::fmax(worst, std::fabs(todbl(abs(got[e] - ref[e]))));
        else if (std::memcmp(&got[e], &ref[e], sizeof(T)) != 0)
            worst = INFINITY;
    }
    if (!(worst <= bound)) {
        failures++;
        std::printf("FAIL %s: %c%c m=%ld n=%ld k=%ld error %g > bound %g\n", what, p.ta, p.tb, p.m, p.n, p.k, worst, bound);
    }
}

int main() {
#ifndef MPLAPACK_GEMM_CUDA_HOST_EMULATION
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) {
        std::printf("no CUDA device: skipped\n");
        return 77;
    }
#endif
    const long sizes[][3] = {{1, 1, 1}, {2, 3, 1}, {7, 5, 3}, {17, 9, 33}, {33, 40, 70}, {16, 16, 16}, {3, 50, 2}};
    const char tr[2] = {'n', 't'};

    for (auto &z : sizes)
        for (char ta : tr)
            for (char tb : tr)
                for (int bm = 0; bm < 3; bm++) {
                    reset_env();
                    const problem p = make(ta, tb, z[0], z[1], z[2], bm, false);
                    const std::vector<T> ref = cpu(p);
                    expect_identical("tiled", p, gpu(p, "tiled"), ref);
                    env("_CUDA_KERNEL", "naive");
                    expect_identical("naive", p, gpu(p, "naive"), ref);
                    if (z[2] > 0) {
                        reset_env();
                        env("_GEMM_OZAKI", "1");
                        const std::vector<T> oref = cpu(p);
                        env("_GEMM_OZAKI", NULL);
                        env("_CUDA_OZAKI", "1");
                        expect_identical("ozaki", p, gpu(p, "ozaki"), oref);
                    }
                }

    const long wsizes[][3] = {{2, 2, 2}, {7, 5, 3}, {17, 9, 33}, {33, 40, 70}};
    const char *cutoffs[] = {"1", "2", "3", "8"};
    for (const char *cut : cutoffs)
        for (auto &z : wsizes) {
            if (cut[0] < '3' && z[0] * z[1] * z[2] > 20000)
                continue;
            for (char ta : tr)
                for (char tb : tr)
                    for (int bm = 0; bm < 3; bm++)
                        for (int integers = 0; integers < 2; integers++) {
                            reset_env();
                            const problem p = make(ta, tb, z[0], z[1], z[2], bm, integers != 0);
                            const std::vector<T> ref = cpu(p);
                            env("_CUDA_WINOGRAD_CUTOFF", cut);
                            const std::vector<T> got = gpu(p, "winograd");
                            if (integers)
                                expect_identical("winograd", p, got, ref);
                            else
                                expect_close("winograd", p, got, ref);
                        }
        }

    // disabled: the call is left to the CPU
    {
        reset_env();
        env("_CUDA", "0");
        problem p = make('n', 'n', 8, 8, 8, 2, false);
        checks++;
        if (Rgemm_gpu(true, true, p.m, p.n, p.k, p.alpha, p.A.data(), p.lda, p.B.data(), p.ldb, p.beta, p.C.data(), p.ldc)) {
            failures++;
            std::printf("FAIL %s_CUDA=0 still ran on the GPU\n", PREFIX);
        }
    }

    std::printf("%s CUDA Rgemm: %ld checks, %ld failures\n", PREFIX, checks, failures);
    return failures == 0 ? 0 : 1;
}

#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
// the GPU code of the library, run on the host
#if defined(RGEMM_TEST_DD)
#include "../../dd/cuda/Rgemm_gpu_cuda.cu"
#else
#include "../../qd/cuda/Rgemm_gpu_cuda.cu"
#endif
#endif
