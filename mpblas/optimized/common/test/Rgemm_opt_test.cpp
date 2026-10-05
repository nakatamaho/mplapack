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
 * Test of the optimized CPU Rgemm of the dd, qd and binary128 libraries
 * (compile with -DRGEMM_TEST_DD, -DRGEMM_TEST_QD or -DRGEMM_TEST_BINARY128
 * and link with the corresponding libmplapack_<type>_opt):
 *
 *  - the blocked engine (and, for dd, its SIMD kernel) is bit-identical to
 *    the plain OpenMP loops openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp, for all
 *    transpose combinations, odd sizes, leading dimensions larger than
 *    needed, alpha/beta in {general, 0, 1};
 *  - Winograd (cutoffs 1, 2, 3, 8), Ozaki and (qd) the branch-free kernel are
 *    exact on small integer matrices and agree with the conventional result
 *    within 2^(20-p) (k |alpha| max|A| max|B| + |beta| max|C|) on random
 *    real matrices (p: precision in bits).
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
static double todbl(const T &x) { return x.x[0]; }
static T absT(const T &x) { return abs(x); }
#elif defined(RGEMM_TEST_QD)
#include <mpblas_qd.h>
typedef qd_real T;
static const char *PREFIX = "MPLAPACK_QD";
static const int PREC = 212;
static double todbl(const T &x) { return x.x[0]; }
static T absT(const T &x) { return abs(x); }
#elif defined(RGEMM_TEST_BINARY128)
#include <mpblas_binary128.h>
typedef mplapack_binary128_t T;
static const char *PREFIX = "MPLAPACK_BINARY128";
static const int PREC = 113;
static double todbl(const T &x) { return (double)x; }
static T absT(const T &x) { return x < 0 ? -x : x; }
#else
#error "define RGEMM_TEST_DD, RGEMM_TEST_QD or RGEMM_TEST_BINARY128"
#endif

void Rgemm_NN_omp(mplapackint m, mplapackint n, mplapackint k, T alpha, T *A, mplapackint lda, T *B, mplapackint ldb, T beta, T *C, mplapackint ldc);
void Rgemm_NT_omp(mplapackint m, mplapackint n, mplapackint k, T alpha, T *A, mplapackint lda, T *B, mplapackint ldb, T beta, T *C, mplapackint ldc);
void Rgemm_TN_omp(mplapackint m, mplapackint n, mplapackint k, T alpha, T *A, mplapackint lda, T *B, mplapackint ldb, T beta, T *C, mplapackint ldc);
void Rgemm_TT_omp(mplapackint m, mplapackint n, mplapackint k, T alpha, T *A, mplapackint lda, T *B, mplapackint ldb, T beta, T *C, mplapackint ldc);

static std::mt19937_64 rng(20261005);
static double u11() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }

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

static void setenv_s(const char *name, const char *value) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%s%s", PREFIX, name);
    if (value)
        setenv(buf, value, 1);
    else
        unsetenv(buf);
}

static void reset_env() {
    setenv_s("_GEMM_BLOCKED", NULL);
    setenv_s("_GEMM_MIN_MNK", "1");
    setenv_s("_GEMM_WINOGRAD_CUTOFF", NULL);
    setenv_s("_GEMM_OZAKI", NULL);
    setenv_s("_GEMM_OZAKI_SLICES", NULL);
    unsetenv("MPLAPACK_QD_GEMM_BF");
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

static std::vector<T> run_rgemm(const problem &p) {
    problem q = p;
    const char ta[2] = {p.ta, 0}, tb[2] = {p.tb, 0};
    Rgemm(ta, tb, q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc);
    return q.C;
}

static std::vector<T> run_omp(const problem &p) {
    problem q = p;
    if (p.ta == 'n' && p.tb == 'n')
        Rgemm_NN_omp(q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc);
    else if (p.ta == 'n')
        Rgemm_NT_omp(q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc);
    else if (p.tb == 'n')
        Rgemm_TN_omp(q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc);
    else
        Rgemm_TT_omp(q.m, q.n, q.k, q.alpha, q.A.data(), q.lda, q.B.data(), q.ldb, q.beta, q.C.data(), q.ldc);
    return q.C;
}

static double maxabs(const std::vector<T> &v) {
    double r = 0.0;
    for (const auto &x : v)
        r = std::fmax(r, std::fabs(todbl(x)));
    return r;
}

static long failures = 0, checks = 0;

static void expect_identical(const char *what, const problem &p, const std::vector<T> &got, const std::vector<T> &ref) {
    checks++;
    if (std::memcmp(got.data(), ref.data(), sizeof(T) * got.size()) != 0) {
        failures++;
        std::printf("FAIL %s: %c%c m=%ld n=%ld k=%ld not bit-identical\n", what, p.ta, p.tb, p.m, p.n, p.k);
    }
}

static void expect_close(const char *what, const problem &p, const std::vector<T> &got, const std::vector<T> &ref) {
    checks++;
    const double bound = std::ldexp(1.0, 20 - PREC) * (p.k * std::fabs(todbl(p.alpha)) * maxabs(p.A) * maxabs(p.B) + std::fabs(todbl(p.beta)) * maxabs(p.C));
    double worst = 0.0;
    for (long j = 0; j < p.n; j++)
        for (long i = 0; i < p.m; i++) {
            const T d = absT(got[i + j * p.ldc] - ref[i + j * p.ldc]);
            worst = std::fmax(worst, todbl(d));
        }
    // entries outside the m x n block must be untouched
    for (size_t e = 0; e < got.size(); e++) {
        const long i = (long)(e % p.ldc), j = (long)(e / p.ldc);
        if ((i >= p.m || j >= p.n) && std::memcmp(&got[e], &ref[e], sizeof(T)) != 0)
            worst = INFINITY;
    }
    if (!(worst <= bound)) {
        failures++;
        std::printf("FAIL %s: %c%c m=%ld n=%ld k=%ld error %g > bound %g\n", what, p.ta, p.tb, p.m, p.n, p.k, worst, bound);
    }
}

int main() {
    const long sizes[][3] = {{1, 1, 1}, {2, 3, 1}, {7, 5, 3}, {17, 9, 33}, {65, 70, 300}, {130, 33, 257}, {64, 64, 64}, {3, 200, 2}, {41, 37, 0}};
    const char tr[2] = {'n', 't'};

    // Conventional paths: bit-identical to the OpenMP loops.
    reset_env();
    for (auto &z : sizes)
        for (char ta : tr)
            for (char tb : tr)
                for (int bm = 0; bm < 3; bm++) {
                    const problem p = make(ta, tb, z[0], z[1], z[2], bm, false);
                    const std::vector<T> ref = run_omp(p);
                    reset_env();
                    expect_identical("blocked", p, run_rgemm(p), ref);
                    setenv_s("_GEMM_BLOCKED", "0");
                    expect_identical("unblocked", p, run_rgemm(p), ref);
                    reset_env();
                    setenv_s("_GEMM_MIN_MNK", "1000000000");
                    expect_identical("small-size path", p, run_rgemm(p), ref);
                    reset_env();
                }

    // Winograd, Ozaki, branch-free qd: exact on integers, bounded on reals.
    const long wsizes[][3] = {{1, 1, 1}, {7, 5, 3}, {17, 9, 33}, {64, 64, 64}, {65, 70, 129}, {100, 3, 40}};
    for (int variant = 0; variant < 7; variant++) {
        const char *name = "?";
        reset_env();
        switch (variant) {
        case 0: name = "winograd cutoff 1"; setenv_s("_GEMM_WINOGRAD_CUTOFF", "1"); break;
        case 1: name = "winograd cutoff 2"; setenv_s("_GEMM_WINOGRAD_CUTOFF", "2"); break;
        case 2: name = "winograd cutoff 3"; setenv_s("_GEMM_WINOGRAD_CUTOFF", "3"); break;
        case 3: name = "winograd cutoff 8"; setenv_s("_GEMM_WINOGRAD_CUTOFF", "8"); break;
        case 4: name = "ozaki"; setenv_s("_GEMM_OZAKI", "1"); break;
        case 5: name = "ozaki + winograd"; setenv_s("_GEMM_OZAKI", "1"); setenv_s("_GEMM_WINOGRAD_CUTOFF", "4"); break;
        case 6:
#if defined(RGEMM_TEST_QD)
            name = "qd branch-free";
            setenv("MPLAPACK_QD_GEMM_BF", "1", 1);
            break;
#else
            continue;
#endif
        }
        for (auto &z : wsizes) {
            // small cutoffs recurse deeply: keep them to small sizes
            if (variant <= 2 && z[0] * z[1] * z[2] > 20000)
                continue;
            for (char ta : tr)
                for (char tb : tr)
                    for (int bm = 0; bm < 3; bm++)
                        for (int integers = 0; integers < 2; integers++) {
                            const problem p = make(ta, tb, z[0], z[1], z[2], bm, integers != 0);
                            const std::vector<T> ref = run_omp(p);
                            const std::vector<T> got = run_rgemm(p);
                            if (integers)
                                expect_identical(name, p, got, ref);
                            else
                                expect_close(name, p, got, ref);
                        }
        }
    }

    // Ozaki falls back to the conventional result on non-finite input.
    {
        reset_env();
        problem p = make('n', 'n', 20, 20, 20, 2, false);
        p.A[5] = (T)INFINITY;
        const std::vector<T> ref = run_omp(p);
        setenv_s("_GEMM_OZAKI", "1");
        const std::vector<T> got = run_rgemm(p);
        checks++;
        bool same = true;
        for (size_t e = 0; e < got.size(); e++)
            if (std::memcmp(&got[e], &ref[e], sizeof(T)) != 0 && !(std::isnan(todbl(got[e])) && std::isnan(todbl(ref[e]))))
                same = false;
        if (!same) {
            failures++;
            std::printf("FAIL ozaki fallback on non-finite input\n");
        }
    }

    std::printf("%s: %ld checks, %ld failures\n", PREFIX, checks, failures);
    return failures == 0 ? 0 : 1;
}
