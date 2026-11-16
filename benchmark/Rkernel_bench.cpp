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
 * Kernel benchmark: the building blocks of the LU / Cholesky factorizations,
 * solvers and inverses (GEMM, SYRK/HERK, TRSM, the getrf2 / potrf2 panels,
 * laswp, trtri, lauum), timed alone, on the CPU or on the GPU.  One source,
 * three programs (see kernel_bench.py and README.kernel_bench.md):
 *
 *   Rkernel_bench.<b>_opt          CPU: the routines of libmplapack_<b>_opt
 *                                  (OpenMP), wall-clock time of the call.
 *   Rkernel_bench.<b>_cuda_kernel  dd/td/qd GPU (Rkernel_bench_<b>_cuda.cu,
 *                                  KERNEL_BENCH_BF_CUDA): the device kernels
 *                                  of mpblas/optimized/common and
 *                                  mplapack/optimized/common, compiled into
 *                                  this program from the library's internal
 *                                  headers, on matrices already on the
 *                                  device; CUDA events around the kernels.
 *   Rkernel_bench.mpfr_cuda_kernel MPFR GPU (KERNEL_BENCH_MPFR_CUDA): the
 *                                  entry points <R>_mpfr_cuda of
 *                                  libmplapack_mpfr_opt_cuda.  Their kernels
 *                                  cannot be called on device data from
 *                                  outside the library, so the time is the
 *                                  span of the call's kernels recorded by
 *                                  CUPTI (start of the first kernel to the
 *                                  end of the last one): packing into
 *                                  cu_freal / cu_mpfr and the transfers are
 *                                  not included.
 *
 * usage: Rkernel_bench.X KERNEL FORM SIZES [options]
 *   KERNEL FORM    gemm NN|NT|TN, gemm_update NN (k = 64), syrk LN|UT[_k64],
 *                  trsm LLNU|RLTN|LUTN|LUNN[_p64], getrf2 P64, potrf2 L,
 *                  laswp R64, trtri UN, lauum L; complex: cgemm NN|NC|CN,
 *                  cherk LN|UC[_k64], ctrsm LLNU|RLCN|LUCN|LUNN[_p64]
 *   SIZES          comma separated list of n
 *   --label L      type label of the output (default dd, td, qd, mpfr<prec>)
 *   --reps R       best of R measured runs (default 3)
 *   --long S       one measured run when the warm-up takes more than S
 *                  seconds (default 10)
 *   --warmup W     warm-up runs (default 1 on the GPU, 0 on the CPU)
 *   --dump DIR     write the results of the sizes in --dump-sizes to DIR
 *   --dump-sizes   comma separated (default 64,200)
 * Output (CSV, no header): label,kernel,form,n,method,seconds,gflops,runs
 * The MPFR programs need MPFRXX_DEFAULT_PRECISION_BITS (all threads).
 */

#if defined(KERNEL_BENCH_BF_CUDA)
// nvcc already provides the two-argument long min/max overloads.
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#endif
#include <mpblas.h>
#include <mplapack.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#if defined(KERNEL_BENCH_BF_CUDA)
#pragma GCC visibility push(hidden)
#include "../mplapack/optimized/common/Rgetrf_gpu_common.cuh"
#include "../mplapack/optimized/common/Rinverse_gpu_common.cuh"
#pragma GCC visibility pop
#endif

#if defined(KERNEL_BENCH_MPFR_CUDA)
#include <mpblas_mpfr_cuda.h>
#include <cuda_runtime_api.h>
#include <cupti.h>
#endif

namespace {

#if defined(MPLAPACK_BUILD_WITH_DD)
const char *default_label() { return "dd"; }
#elif defined(MPLAPACK_BUILD_WITH_TD)
const char *default_label() { return "td"; }
#elif defined(MPLAPACK_BUILD_WITH_QD)
const char *default_label() { return "qd"; }
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
std::string mpfr_label;
const char *default_label() {
    mpfr_label = "mpfr" + std::to_string((long)mpfr_get_default_prec());
    return mpfr_label.c_str();
}
#else
#error "Rkernel_bench: dd, td, qd or mpfr"
#endif

#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_CUDA)
const char *method_name = "gpu";
#else
const char *method_name = "cpu";
#endif

// ---- deterministic full-precision data ----

uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
// 53-bit uniform value in [-1, 1)
double uniform(uint64_t &state) {
    state = mix(state);
    return (double)(int64_t)(state >> 11) * 0x1p-52 - 1.0;
}

#if defined(MPLAPACK_BUILD_WITH_MPFR)
// MPFR: every bit of the significand is random (mpfr_urandomb), in [-1, 1).
struct generator {
    gmp_randstate_t rs;
    explicit generator(uint64_t seed) {
        gmp_randinit_default(rs);
        gmp_randseed_ui(rs, (unsigned long)seed);
    }
    ~generator() { gmp_randclear(rs); }
    REAL real() {
        REAL x;
        mpfr_urandomb(x.get_mpfr_t(), rs);
        mpfr_mul_2ui(x.get_mpfr_t(), x.get_mpfr_t(), 1, MPFR_RNDN);
        mpfr_sub_ui(x.get_mpfr_t(), x.get_mpfr_t(), 1, MPFR_RNDN);
        return x;
    }
};
#else
// dd/td/qd: each limb random.
struct generator {
    uint64_t state;
    explicit generator(uint64_t seed) : state(mix(seed)) {}
    REAL real() {
        const int limbs = (int)(sizeof(REAL) / sizeof(double));
        REAL x = uniform(state);
        double scale = 1.0;
        for (int q = 1; q < limbs; q++) {
            scale *= 0x1p-53;
            x += REAL(uniform(state) * scale);
        }
        return x;
    }
};
#endif

template <class T> T draw(generator &g);
template <> REAL draw<REAL>(generator &g) { return g.real(); }
template <> COMPLEX draw<COMPLEX>(generator &g) {
    REAL re = g.real();
    REAL im = g.real();
    return COMPLEX(re, im);
}
REAL magnitude(const REAL &x) { return abs(x); }
REAL magnitude(const COMPLEX &x) { return abs(x.real()) + abs(x.imag()); }

// Column-major host matrix.
template <class T> struct matrix {
    long rows = 0, cols = 0;
    std::vector<T> v;
    void resize(long r, long c) {
        rows = r;
        cols = c;
        v.assign((size_t)std::max(1L, r) * (size_t)std::max(1L, c), T());
    }
    T *data() { return v.data(); }
    long ld() const { return std::max(1L, rows); }
    T &at(long i, long j) { return v[(size_t)i + (size_t)j * ld()]; }
};

// kind G: uniform; T: triangular factor, off-diagonal / rows, diagonal in [1, 1.5)
// (well conditioned for unit and non-unit solves and inverses); S: symmetric
// positive definite, diagonal rows (diagonally dominant).
template <class T> void fill(matrix<T> &a, long r, long c, char kind, uint64_t seed) {
    a.resize(r, c);
    generator g(seed);
    for (long j = 0; j < c; j++)
        for (long i = 0; i < r; i++) {
            T x = draw<T>(g);
            if (kind == 'T' && i != j)
                x = x / REAL((double)r);
            if (kind == 'T' && i == j)
                x = T(REAL(1) + magnitude(x) / REAL(2));
            if (kind == 'S' && i == j)
                x = T(REAL((double)r));
            a.at(i, j) = x;
        }
    if (kind == 'S')
        for (long j = 0; j < c; j++)
            for (long i = 0; i < j && i < r; i++)
                a.at(i, j) = a.at(j, i);
}

// leading r x c block of src -> dst (same leading dimensions)
template <class T> void restore(matrix<T> &dst, matrix<T> &src, long r, long c) {
    for (long j = 0; j < c; j++)
        for (long i = 0; i < r; i++)
            dst.at(i, j) = src.at(i, j);
}

// ---- hexadecimal dumps (bit comparison) ----

void hex_limbs(FILE *f, const REAL &x) {
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    mpfr_fprintf(f, "%Ra", x.get_mpfr_t());
#else
    const int limbs = (int)(sizeof(REAL) / sizeof(double));
    for (int q = 0; q < limbs; q++)
        std::fprintf(f, q ? ",%a" : "%a", x.x[q]);
#endif
}
void hex(FILE *f, const REAL &x) { hex_limbs(f, x); }
void hex(FILE *f, const COMPLEX &x) {
    hex_limbs(f, x.real());
    std::fputc(';', f);
    hex_limbs(f, x.imag());
}
template <class T> void dump(const std::string &file, matrix<T> &a, long r, long c, const std::vector<mplapackint> *piv, long info) {
    FILE *f = std::fopen(file.c_str(), "w");
    if (!f)
        throw std::runtime_error("cannot write " + file);
    std::fprintf(f, "info %ld\n", info);
    if (piv)
        for (size_t i = 0; i < piv->size(); i++)
            std::fprintf(f, "ipiv %ld\n", (long)(*piv)[i]);
    for (long j = 0; j < c; j++)
        for (long i = 0; i < r; i++) {
            hex(f, a.at(i, j));
            std::fputc('\n', f);
        }
    std::fclose(f);
}

// ---- flop counts (LAPACK Working Note 41 / the usual BLAS counts) ----

double flops_getrf(double m, double n) {
    double muls, adds;
    if (m < n) {
        muls = 0.5 * m * (m * (n - m / 3 - 1) + n) + 2 * m / 3;
        adds = 0.5 * m * (m * (n - m / 3) - n) + m / 6;
    } else {
        muls = 0.5 * n * (n * (m - n / 3 - 1) + m) + 2 * n / 3;
        adds = 0.5 * n * (n * (m - n / 3) - m) + n / 6;
    }
    return muls + adds;
}
double flops_potrf(double n) { return n * n * n / 3 + n * n / 2 + n / 6; }
double flops_trtri(double n) { return n * n * n / 3 + 2 * n / 3; }

// ---- the case being timed ----

struct spec {
    std::string kernel, form;
    bool complex = false;
    // derived
    char side = 'L', uplo = 'L', trans = 'N', diag = 'N', ta = 'N', tb = 'N';
    bool p64 = false, k64 = false;
};

spec parse(const std::string &kernel, const std::string &form) {
    spec s;
    s.kernel = kernel;
    s.form = form;
    s.complex = kernel[0] == 'c';
    std::string f = form;
    const size_t us = f.find('_');
    if (us != std::string::npos) {
        const std::string suffix = f.substr(us + 1);
        f = f.substr(0, us);
        if (suffix == "p64")
            s.p64 = true;
        else if (suffix == "k64")
            s.k64 = true;
        else
            throw std::runtime_error("unknown form suffix " + suffix);
    }
    const std::string k = s.complex ? kernel.substr(1) : kernel;
    if (k == "gemm" || k == "gemm_update") {
        if (f.size() != 2)
            throw std::runtime_error("gemm form XY");
        s.ta = f[0], s.tb = f[1];
    } else if (k == "syrk" || k == "herk") {
        if (f.size() != 2)
            throw std::runtime_error("syrk form UPLO TRANS");
        s.uplo = f[0], s.trans = f[1];
    } else if (k == "trsm") {
        if (f.size() != 4)
            throw std::runtime_error("trsm form SIDE UPLO TRANS DIAG");
        s.side = f[0], s.uplo = f[1], s.trans = f[2], s.diag = f[3];
    } else if (k == "trtri") {
        s.uplo = f[0], s.diag = f[1];
    } else if (k == "lauum" || k == "potrf2") {
        s.uplo = f[0];
    } else if (k != "getrf2" && k != "laswp") {
        throw std::runtime_error("unknown kernel " + kernel);
    }
    return s;
}

const long NB = 64; // block width of the factorizations (iMlaenv)

// Operands of one size n: a (input), b (input or in/out), c (in/out).
// The pools have the largest size; a size uses their leading blocks.
template <class T> struct operands {
    matrix<T> a, b, c, c0; // c0: pristine in/out operand
    std::vector<mplapackint> piv;
    long am = 0, an = 0, bm = 0, bn = 0, cm = 0, cn = 0; // used blocks at size n
    bool has_b = false, b_out = false;
};

double flops(const spec &s, long n) {
    const double x = (double)n, f = s.complex ? 4.0 : 1.0;
    const std::string k = s.complex ? s.kernel.substr(1) : s.kernel;
    if (k == "gemm")
        return 2.0 * f * x * x * x;
    if (k == "gemm_update")
        return 2.0 * f * x * x * NB;
    if (k == "syrk" || k == "herk") {
        const double kk = s.k64 ? NB : x;
        return f * kk * x * (x + 1);
    }
    if (k == "trsm") {
        // B is m x n; Left: n m^2, Right: m n^2
        const double tri = s.p64 ? NB : x, other = x;
        return f * other * tri * tri;
    }
    if (k == "getrf2")
        return flops_getrf(x, NB);
    if (k == "potrf2" || k == "lauum")
        return flops_potrf(x);
    if (k == "trtri")
        return flops_trtri(x);
    return 0.0; // laswp
}

// ---- timing of one call ----

#if defined(KERNEL_BENCH_MPFR_CUDA)
// CUPTI: span of the kernels of a call.
uint64_t kernel_first = UINT64_MAX, kernel_last = 0;
long kernel_count = 0;
void CUPTIAPI buffer_requested(uint8_t **buffer, size_t *size, size_t *max_records) {
    *size = 8 << 20;
    *buffer = static_cast<uint8_t *>(std::aligned_alloc(8, *size));
    *max_records = 0;
}
void CUPTIAPI buffer_completed(CUcontext, uint32_t, uint8_t *buffer, size_t, size_t valid) {
    CUpti_Activity *record = NULL;
    while (cuptiActivityGetNextRecord(buffer, valid, &record) == CUPTI_SUCCESS) {
        if (record->kind == CUPTI_ACTIVITY_KIND_CONCURRENT_KERNEL || record->kind == CUPTI_ACTIVITY_KIND_KERNEL) {
            const CUpti_ActivityKernel11 *k = reinterpret_cast<const CUpti_ActivityKernel11 *>(record);
            kernel_first = std::min<uint64_t>(kernel_first, k->start);
            kernel_last = std::max<uint64_t>(kernel_last, k->end);
            kernel_count++;
        }
    }
    std::free(buffer);
}
void cupti_check(CUptiResult r, const char *what) {
    if (r != CUPTI_SUCCESS) {
        const char *msg = "";
        cuptiGetResultString(r, &msg);
        throw std::runtime_error(std::string("CUPTI ") + what + ": " + msg);
    }
}
void cupti_start() {
    cupti_check(cuptiActivityRegisterCallbacks(buffer_requested, buffer_completed), "register");
    cupti_check(cuptiActivityEnable(CUPTI_ACTIVITY_KIND_CONCURRENT_KERNEL), "enable");
}
template <class Call> double timed(Call call) {
    cudaDeviceSynchronize();
    cupti_check(cuptiActivityFlushAll(1), "flush");
    kernel_first = UINT64_MAX, kernel_last = 0, kernel_count = 0;
    if (!call())
        throw std::runtime_error("entry point returned false (not computed on the GPU)");
    if (cudaDeviceSynchronize() != cudaSuccess)
        throw std::runtime_error("CUDA error after the call");
    cupti_check(cuptiActivityFlushAll(1), "flush");
    if (kernel_count == 0)
        throw std::runtime_error("no kernel recorded");
    return (double)(kernel_last - kernel_first) * 1e-9;
}
#elif defined(KERNEL_BENCH_BF_CUDA)
void cuda_check(cudaError_t e, const char *what) {
    if (e != cudaSuccess)
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(e));
}
// The kernels synchronize before returning (blas_device, gemm_device,
// job_device, the stream contexts' finish()).
template <class Call> double timed(Call call) {
    cudaEvent_t start, stop;
    cuda_check(cudaEventCreate(&start), "event");
    cuda_check(cudaEventCreate(&stop), "event");
    cuda_check(cudaDeviceSynchronize(), "synchronize");
    cuda_check(cudaEventRecord(start, 0), "record");
    cuda_check(cudaEventSynchronize(start), "synchronize start");
    const int e = call();
    cuda_check(cudaEventRecord(stop, 0), "record");
    cuda_check(cudaEventSynchronize(stop), "synchronize stop");
    if (e)
        throw std::runtime_error("kernel error " + std::to_string(e));
    float ms = 0;
    cuda_check(cudaEventElapsedTime(&ms, start, stop), "elapsed");
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return ms * 1e-3;
}
#else
template <class Call> double timed(Call call) {
    auto t0 = std::chrono::steady_clock::now();
    call();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
#endif

// ---- dd/td/qd GPU: device copies and the kernels ----

#if defined(KERNEL_BENCH_BF_CUDA)
#if defined(MPLAPACK_BUILD_WITH_DD)
using Ops = mplapack_gemm::dd_ops;
#elif defined(MPLAPACK_BUILD_WITH_TD)
using Ops = mplapack_gemm::td_ops;
#else
using Ops = mplapack_gemm::qd_ops;
#endif
using F = Ops::F;
using DR = mplapack_blas_dev::dreal<Ops>;
using DC = mplapack_blas_dev::dcomplex<DR>;
const int LIMBS = (int)(sizeof(REAL) / sizeof(double));
static_assert(sizeof(F) == sizeof(REAL), "device and host limbs");
static_assert(sizeof(DR) == sizeof(F) && sizeof(DC) == 2 * sizeof(F), "device element layout");

void pack(const REAL &x, double *d) {
    for (int q = 0; q < LIMBS; q++)
        d[q] = x.x[q];
}
void pack(const COMPLEX &x, double *d) {
    pack(x.real(), d);
    pack(x.imag(), d + LIMBS);
}
void unpack(const double *d, REAL &x) {
    for (int q = 0; q < LIMBS; q++)
        x.x[q] = d[q];
}
void unpack(const double *d, COMPLEX &x) {
    REAL re, im;
    unpack(d, re);
    unpack(d + LIMBS, im);
    x = COMPLEX(re, im);
}

// Device copy of the leading r x c block of a host matrix, ld = r.
template <class T> struct device_matrix {
    static const int width = (int)(sizeof(T) / sizeof(REAL)) * LIMBS; // doubles per element
    double *d = nullptr;
    long r = 0, c = 0;
    std::vector<double> staging;
    ~device_matrix() { cudaFree(d); }
    void allocate(long rows, long cols) {
        cudaFree(d);
        d = nullptr;
        r = rows, c = cols;
        cuda_check(cudaMalloc(&d, (size_t)std::max(1L, r) * std::max(1L, c) * width * sizeof(double)), "cudaMalloc");
    }
    void stage(matrix<T> &h) {
        staging.resize((size_t)r * c * width);
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                pack(h.at(i, j), &staging[((size_t)i + (size_t)j * r) * width]);
    }
    void upload() { cuda_check(cudaMemcpy(d, staging.data(), staging.size() * sizeof(double), cudaMemcpyHostToDevice), "upload"); }
    void download(matrix<T> &h) {
        std::vector<double> t((size_t)r * c * width);
        cuda_check(cudaMemcpy(t.data(), d, t.size() * sizeof(double), cudaMemcpyDeviceToHost), "download");
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                unpack(&t[((size_t)i + (size_t)j * r) * width], h.at(i, j));
    }
    long ld() const { return std::max(1L, r); }
    DR *real() const { return reinterpret_cast<DR *>(d); }
    DC *cplx() const { return reinterpret_cast<DC *>(d); }
};

// Constants on the device: 1 and -1 (as the factorizations pass them).
struct device_constants {
    F *p = nullptr;
    device_constants() {
        F c[2] = {Ops::from_d(1.0), Ops::from_d(-1.0)};
        cuda_check(cudaMalloc(&p, sizeof(c)), "cudaMalloc");
        cuda_check(cudaMemcpy(p, c, sizeof(c), cudaMemcpyHostToDevice), "upload");
    }
    ~device_constants() { cudaFree(p); }
};
#endif

// ---- one kernel at one size ----

template <class T> struct runner {
    spec s;
    operands<T> op;
    long nmax = 0;
    std::string label;
    int reps = 3, warmup = 0;
    double long_run = 10.0;
    std::string dump_dir;
    std::vector<long> dump_sizes;
#if defined(KERNEL_BENCH_BF_CUDA)
    device_constants constants;
#endif

    std::string base() const { return s.complex ? s.kernel.substr(1) : s.kernel; }

    // pools at the largest size
    void prepare() {
        const std::string k = base();
        const uint64_t seed = 20261104;
        const long n = nmax;
        op.has_b = false;
        if (k == "gemm") {
            fill(op.a, n, n, 'G', seed + 1);
            fill(op.b, n, n, 'G', seed + 2);
            op.has_b = true;
            fill(op.c0, n, n, 'G', seed + 3);
        } else if (k == "gemm_update") {
            fill(op.a, n, NB, 'G', seed + 1);
            fill(op.b, NB, n, 'G', seed + 2);
            op.has_b = true;
            fill(op.c0, n, n, 'G', seed + 3);
        } else if (k == "syrk" || k == "herk") {
            const long kk = s.k64 ? NB : n;
            if (s.trans == 'N')
                fill(op.a, n, kk, 'G', seed + 1);
            else
                fill(op.a, kk, n, 'G', seed + 1);
            fill(op.c0, n, n, 'G', seed + 3);
        } else if (k == "trsm") {
            const long tri = s.p64 ? NB : n;
            fill(op.a, tri, tri, 'T', seed + 1);
            if (s.side == 'L')
                fill(op.c0, tri, n, 'G', seed + 3);
            else
                fill(op.c0, n, tri, 'G', seed + 3);
        } else if (k == "getrf2") {
            fill(op.c0, n, NB, 'G', seed + 3);
        } else if (k == "potrf2") {
            fill(op.c0, n, n, 'S', seed + 3);
        } else if (k == "laswp") {
            fill(op.c0, n, n, 'G', seed + 3);
        } else if (k == "trtri" || k == "lauum") {
            fill(op.c0, n, n, 'T', seed + 3);
        }
        op.c = op.c0;
    }

    // the blocks of size n
    void shape(long n) {
        const std::string k = base();
        op.am = op.an = op.bm = op.bn = 0;
        if (k == "gemm")
            op.am = op.an = op.bm = op.bn = op.cm = op.cn = n;
        else if (k == "gemm_update") {
            op.am = n, op.an = NB, op.bm = NB, op.bn = n, op.cm = op.cn = n;
        } else if (k == "syrk" || k == "herk") {
            const long kk = s.k64 ? NB : n;
            op.am = s.trans == 'N' ? n : kk;
            op.an = s.trans == 'N' ? kk : n;
            op.cm = op.cn = n;
        } else if (k == "trsm") {
            const long tri = s.p64 ? NB : n;
            op.am = op.an = tri;
            op.cm = s.side == 'L' ? tri : n;
            op.cn = s.side == 'L' ? n : tri;
        } else if (k == "getrf2") {
            op.cm = n, op.cn = NB;
        } else {
            op.cm = op.cn = n;
        }
        op.piv.assign(0, 0);
        if (k == "getrf2")
            op.piv.assign(std::min(n, NB), 0);
        if (k == "laswp") {
            // 64 interchanges, row i with a row p >= i (as after a panel)
            uint64_t st = 777 + (uint64_t)n;
            const long rows = std::min(n, NB);
            op.piv.resize(rows);
            for (long i = 0; i < rows; i++)
                op.piv[i] = i + 1 + (long)(mix(st++) % (uint64_t)(n - i));
        }
    }

#if defined(KERNEL_BENCH_BF_CUDA)
    device_matrix<T> da, db, dc, dwork;
    long *dpiv = nullptr;
    // contexts (streams, workspace) are created outside the timed region
    std::unique_ptr<mplapack_lapack_gpu::lu_context<mplapack_lapack_gpu::real_lu<Ops>, DR>> lu;
    std::unique_ptr<mplapack_lapack_gpu::real_cholesky<Ops>> chol;
    std::unique_ptr<mplapack_inverse_gpu::context<Ops, false>> inv;
    mplapack_blas_dev::job_Cgemm<Ops> *jg = nullptr;
    mplapack_blas_dev::job_Cherk<Ops> *jh = nullptr;
    mplapack_blas_dev::job_Ctrsm<Ops> *jt = nullptr;
    void device_setup(long n) {
        const std::string k = base();
        if (op.am) {
            da.allocate(op.am, op.an);
            da.stage(op.a);
            da.upload();
        }
        if (op.has_b) {
            db.allocate(op.bm, op.bn);
            db.stage(op.b);
            db.upload();
        }
        dc.allocate(op.cm, op.cn);
        dc.stage(op.c0);
        if (k == "gemm" || k == "gemm_update")
            dwork.allocate(std::max(n, NB), n); // alpha*op(B), k x n
        lu.reset();
        chol.reset();
        inv.reset();
        if (k == "getrf2" && !s.complex) {
            DR sfmin;
            pack(Rlamch("S"), reinterpret_cast<double *>(&sfmin));
            lu.reset(new mplapack_lapack_gpu::lu_context<mplapack_lapack_gpu::real_lu<Ops>, DR>(std::max(op.cm, op.cn), NB, sfmin));
            cuda_check((cudaError_t)lu->error, "LU context");
        }
        if (k == "potrf2" && !s.complex) {
            chol.reset(new mplapack_lapack_gpu::real_cholesky<Ops>(n, NB));
            cuda_check((cudaError_t)chol->error, "Cholesky context");
        }
        if ((k == "trtri" || k == "lauum") && !s.complex) {
            inv.reset(new mplapack_inverse_gpu::context<Ops, false>());
            cuda_check((cudaError_t)inv->error, "inverse context");
        }
        cudaFree(dpiv);
        dpiv = nullptr;
        if (!op.piv.empty()) {
            cuda_check(cudaMalloc(&dpiv, op.piv.size() * sizeof(long)), "cudaMalloc");
            std::vector<long> p(op.piv.begin(), op.piv.end());
            cuda_check(cudaMemcpy(dpiv, p.data(), p.size() * sizeof(long), cudaMemcpyHostToDevice), "upload");
        }
        if (!jg) {
            cuda_check(cudaMalloc(&jg, sizeof(*jg)), "cudaMalloc");
            cuda_check(cudaMalloc(&jh, sizeof(*jh)), "cudaMalloc");
            cuda_check(cudaMalloc(&jt, sizeof(*jt)), "cudaMalloc");
        }
    }
    void before_run() {
        dc.upload();
        if (base() == "laswp" || base() == "getrf2") {
            std::vector<long> p(op.piv.begin(), op.piv.end());
            if (!p.empty())
                cuda_check(cudaMemcpy(dpiv, p.data(), p.size() * sizeof(long), cudaMemcpyHostToDevice), "upload");
        }
    }
    void after_runs() {
        dc.download(op.c);
        if (base() == "getrf2") {
            std::vector<long> p(op.piv.size());
            cuda_check(cudaMemcpy(p.data(), dpiv, p.size() * sizeof(long), cudaMemcpyDeviceToHost), "download");
            for (size_t i = 0; i < p.size(); i++)
                op.piv[i] = p[i];
        }
    }
    long info = 0;
    int call(long n) {
        using namespace mplapack_lapack_gpu;
        const std::string k = base();
        info = 0;
        if (!s.complex) {
            if (k == "gemm" || k == "gemm_update") {
                mplapack_gemm::shape sh{};
                sh.nota = s.ta == 'N', sh.notb = s.tb == 'N';
                sh.m = op.cm, sh.n = op.cn, sh.k = k == "gemm" ? n : NB;
                sh.lda = da.ld(), sh.ldb = db.ld(), sh.ldc = dc.ld();
                mplapack_gemm::gpu_scalars<F> sc{Ops::from_d(-1.0), Ops::from_d(1.0), false, true};
                return mplapack_gemm::gemm_device<Ops>(sh, sc, reinterpret_cast<F *>(da.d), reinterpret_cast<F *>(db.d), reinterpret_cast<F *>(dc.d), reinterpret_cast<F *>(dwork.d));
            }
            mplapack_blas::gpu_lane_args<Ops> g{};
            if (k == "syrk") {
                g.s.op = mplapack_blas::L3_SYRK;
                g.s.upper = s.uplo == 'U', g.s.trans = s.trans != 'N';
                g.s.n = n, g.s.k = s.k64 ? NB : n, g.s.beta_one = 1;
                g.a = reinterpret_cast<F *>(da.d), g.c = reinterpret_cast<F *>(dc.d);
                g.al = constants.p + 1, g.be = constants.p;
                g.lda = da.ld(), g.ldc = dc.ld();
                return mplapack_blas::blas_device<Ops>(g);
            }
            if (k == "trsm") {
                g.s.op = mplapack_blas::L3_TRSM;
                g.s.left = s.side == 'L', g.s.upper = s.uplo == 'U';
                g.s.trans = s.trans != 'N', g.s.nounit = s.diag == 'N';
                g.s.m = op.cm, g.s.n = op.cn, g.s.alpha_one = 1;
                g.a = reinterpret_cast<F *>(da.d), g.c = reinterpret_cast<F *>(dc.d);
                g.al = constants.p, g.lda = da.ld(), g.ldc = dc.ld();
                return mplapack_blas::blas_device<Ops>(g);
            }
            if (k == "getrf2") { // the panel of the device LU (lu() calls recursive())
                info = lu->recursive(op.cm, op.cn, dc.real(), dc.ld(), dpiv);
                lu->finish();
                return lu->error;
            }
            if (k == "potrf2") { // the diagonal block of the device Cholesky
                info = chol_recursive(*chol, s.uplo == 'U', n, dc.real(), dc.ld());
                chol->finish();
                return chol->error;
            }
            if (k == "laswp")
                return laswp<DR>(n, dc.real(), dc.ld(), 1, (long)op.piv.size(), dpiv, 1);
            if (k == "trtri" || k == "lauum") { // as Rinverse_cuda_entry.cuh
                if (k == "trtri") {
                    const char opts[3] = {s.uplo, s.diag, 0};
                    info = inv->tri(s.uplo == 'U', s.diag == 'U', n, iMlaenv(1, "Rtrtri", opts, n, -1, -1, -1), dc.real(), dc.ld());
                } else {
                    const char opts[2] = {s.uplo, 0};
                    inv->lauum(s.uplo == 'U', n, iMlaenv(1, "Rlauum", opts, n, -1, -1, -1), dc.real(), dc.ld());
                }
                inv->finish();
                return inv->error;
            }
        } else {
            using namespace mplapack_blas_dev;
            if (k == "gemm") {
                job_Cgemm<Ops> j{};
                j.ta[0] = s.ta, j.tb[0] = s.tb;
                j.m = op.cm, j.n = op.cn, j.k = n;
                j.lda = da.ld(), j.ldb = db.ld(), j.ldc = dc.ld();
                j.alpha = DC(-1.0, 0.0), j.beta = DC(1.0, 0.0);
                j.A = da.cplx(), j.B = db.cplx(), j.C = dc.cplx();
                return job_device(j, j.units(), jg);
            }
            if (k == "herk") {
                job_Cherk<Ops> j{};
                j.uplo[0] = s.uplo, j.tr[0] = s.trans;
                j.n = n, j.k = s.k64 ? NB : n, j.lda = da.ld(), j.ldc = dc.ld();
                j.alpha = DR(-1.0), j.beta = DR(1.0);
                j.A = da.cplx(), j.C = dc.cplx();
                return job_device(j, j.units(), jh);
            }
            if (k == "trsm") {
                job_Ctrsm<Ops> j{};
                j.side[0] = s.side, j.uplo[0] = s.uplo, j.ta[0] = s.trans, j.diag[0] = s.diag;
                j.m = op.cm, j.n = op.cn, j.lda = da.ld(), j.ldb = dc.ld();
                j.alpha = DC(1.0, 0.0);
                j.A = da.cplx(), j.B = dc.cplx();
                return job_device(j, j.units(), jt);
            }
        }
        throw std::runtime_error("no GPU kernel for " + s.kernel);
    }
#else
    long info = 0;
    void device_setup(long) {}
    void before_run() { restore(op.c, op.c0, op.cm, op.cn); }
    void after_runs() {}
#if defined(KERNEL_BENCH_MPFR_CUDA)
#define KB_CALL(R, ...) return R##_mpfr_cuda(__VA_ARGS__)
#define KB_FACTOR(R, ...) return R##_mpfr_cuda(__VA_ARGS__)
#else
#define KB_CALL(R, ...)                                                                                                                              \
    do {                                                                                                                                             \
        R(__VA_ARGS__);                                                                                                                              \
        return true;                                                                                                                                 \
    } while (0)
#define KB_FACTOR(R, ...) KB_CALL(R, __VA_ARGS__)
#endif
    bool call(long n) {
        const std::string k = base();
        const char sd[2] = {s.side, 0}, up[2] = {s.uplo, 0}, tr[2] = {s.trans, 0}, dg[2] = {s.diag, 0};
        const char ta[2] = {s.ta, 0}, tb[2] = {s.tb, 0};
        info = 0;
        mplapackint inf = 0;
        if constexpr (!std::is_same<T, COMPLEX>::value) {
            if (k == "gemm" || k == "gemm_update")
                KB_CALL(Rgemm, ta, tb, op.cm, op.cn, k == "gemm" ? n : NB, REAL(-1), op.a.data(), op.a.ld(), op.b.data(), op.b.ld(), REAL(1), op.c.data(), op.c.ld());
            if (k == "syrk")
                KB_CALL(Rsyrk, up, tr, n, s.k64 ? NB : n, REAL(-1), op.a.data(), op.a.ld(), REAL(1), op.c.data(), op.c.ld());
            if (k == "trsm")
                KB_CALL(Rtrsm, sd, up, tr, dg, op.cm, op.cn, REAL(1), op.a.data(), op.a.ld(), op.c.data(), op.c.ld());
            if (k == "getrf2") {
#if defined(KERNEL_BENCH_MPFR_CUDA)
                // nb = iMlaenv = 64 >= min(m, 64): the entry runs the recursive panel
                const bool ok = Rgetrf_mpfr_cuda(op.cm, op.cn, op.c.data(), op.c.ld(), op.piv.data(), inf);
#else
                Rgetrf2(op.cm, op.cn, op.c.data(), op.c.ld(), op.piv.data(), inf);
                const bool ok = true;
#endif
                info = inf;
                return ok;
            }
            if (k == "potrf2") {
#if defined(KERNEL_BENCH_MPFR_CUDA)
                // n <= nb: recursive diagonal block; n > nb: blocked Rpotrf
                const bool ok = Rpotrf_mpfr_cuda(up, n, op.c.data(), op.c.ld(), inf);
#else
                Rpotrf2(up, n, op.c.data(), op.c.ld(), inf);
                const bool ok = true;
#endif
                info = inf;
                return ok;
            }
#if !defined(KERNEL_BENCH_MPFR_CUDA)
            if (k == "laswp") {
                Rlaswp(n, op.c.data(), op.c.ld(), 1, (mplapackint)op.piv.size(), op.piv.data(), 1);
                return true;
            }
            if (k == "trtri") {
                Rtrtri(up, dg, n, op.c.data(), op.c.ld(), inf);
                info = inf;
                return true;
            }
            if (k == "lauum") {
                Rlauum(up, n, op.c.data(), op.c.ld(), inf);
                info = inf;
                return true;
            }
#endif
        } else {
            if (k == "gemm")
                KB_CALL(Cgemm, ta, tb, op.cm, op.cn, n, COMPLEX(-1, 0), op.a.data(), op.a.ld(), op.b.data(), op.b.ld(), COMPLEX(1, 0), op.c.data(), op.c.ld());
            if (k == "herk")
                KB_CALL(Cherk, up, tr, n, s.k64 ? NB : n, REAL(-1), op.a.data(), op.a.ld(), REAL(1), op.c.data(), op.c.ld());
            if (k == "trsm")
                KB_CALL(Ctrsm, sd, up, tr, dg, op.cm, op.cn, COMPLEX(1, 0), op.a.data(), op.a.ld(), op.c.data(), op.c.ld());
        }
        throw std::runtime_error("no " + std::string(method_name) + " routine for " + s.kernel + " " + s.form);
    }
#undef KB_CALL
#undef KB_FACTOR
#endif

    void run(long n) {
        shape(n);
        device_setup(n);
        double best = 1e300;
        int runs = 0, nrep = reps;
        for (int w = 0; w < warmup; w++) {
            before_run();
            if (timed([&]() { return call(n); }) > long_run) { // one measured run is enough
                nrep = 1;
                break;
            }
        }
        for (int r = 0; r < nrep; r++) {
            before_run();
            const double t = timed([&]() { return call(n); });
            best = std::min(best, t);
            runs++;
            if (t > long_run)
                break;
        }
        after_runs();
        if (info)
            throw std::runtime_error("info = " + std::to_string(info) + " at n = " + std::to_string(n));
        const double fl = flops(s, n);
        if (fl > 0)
            std::printf("%s,%s,%s,%ld,%s,%.6e,%.6g,%d\n", label.c_str(), s.kernel.c_str(), s.form.c_str(), n, method_name, best, fl / best * 1e-9, runs);
        else
            std::printf("%s,%s,%s,%ld,%s,%.6e,,%d\n", label.c_str(), s.kernel.c_str(), s.form.c_str(), n, method_name, best, runs);
        std::fflush(stdout);
        if (!dump_dir.empty() && std::find(dump_sizes.begin(), dump_sizes.end(), n) != dump_sizes.end()) {
            const std::string file = dump_dir + "/" + label + "_" + s.kernel + "_" + s.form + "_" + std::to_string(n) + "_" + method_name + ".txt";
            dump(file, op.c, op.cm, op.cn, op.piv.empty() || base() == "laswp" ? nullptr : &op.piv, info);
        }
    }
};

std::vector<long> list(const char *text) {
    std::vector<long> v;
    for (const char *p = text; *p;) {
        char *end = nullptr;
        const long x = std::strtol(p, &end, 10);
        if (end == p)
            throw std::runtime_error(std::string("bad list ") + text);
        v.push_back(x);
        p = *end == ',' ? end + 1 : end;
    }
    return v;
}

template <class T> int bench(const spec &s, const std::vector<long> &sizes, int argc, char **argv) {
    runner<T> r;
    r.s = s;
    r.label = default_label();
    r.dump_sizes = {64, 200};
#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_CUDA)
    r.warmup = 1;
#endif
    for (int i = 4; i < argc; i++) {
        const std::string a = argv[i];
        const char *v = i + 1 < argc ? argv[i + 1] : "";
        if (a == "--label")
            r.label = v, i++;
        else if (a == "--reps")
            r.reps = std::atoi(v), i++;
        else if (a == "--long")
            r.long_run = std::atof(v), i++;
        else if (a == "--warmup")
            r.warmup = std::atoi(v), i++;
        else if (a == "--dump")
            r.dump_dir = v, i++;
        else if (a == "--dump-sizes")
            r.dump_sizes = list(v), i++;
        else
            throw std::runtime_error("unknown option " + a);
    }
    r.nmax = *std::max_element(sizes.begin(), sizes.end());
    auto t0 = std::chrono::steady_clock::now();
    r.prepare();
    std::fprintf(stderr, "# %s %s %s: data for n <= %ld in %.1f s\n", r.label.c_str(), s.kernel.c_str(), s.form.c_str(), r.nmax,
                 std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    for (long n : sizes) {
        try {
            r.run(n);
        } catch (const std::exception &e) {
            std::printf("%s,%s,%s,%ld,%s,,,0\n", r.label.c_str(), s.kernel.c_str(), s.form.c_str(), n, method_name);
            std::fflush(stdout);
            std::fprintf(stderr, "# %s %s %s n=%ld FAILED: %s\n", r.label.c_str(), s.kernel.c_str(), s.form.c_str(), n, e.what());
        }
    }
    return 0;
}

} // namespace

#include "Rkernel_new.h"

int main(int argc, char **argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: %s KERNEL FORM SIZES [--label L] [--reps R] [--long S] [--warmup W] [--dump DIR] [--dump-sizes N,N]\n", argv[0]);
        return 2;
    }
    try {
#if defined(MPLAPACK_BUILD_WITH_MPFR)
        if (!std::getenv("MPFRXX_DEFAULT_PRECISION_BITS"))
            throw std::runtime_error("set MPFRXX_DEFAULT_PRECISION_BITS (the precision of every thread)");
        {
            REAL first; // gmpfrxx_mkII reads the environment at the first value of a thread
            if ((long)mpfr_get_default_prec() != std::atol(std::getenv("MPFRXX_DEFAULT_PRECISION_BITS")))
                throw std::runtime_error("the MPFR default precision is not MPFRXX_DEFAULT_PRECISION_BITS");
        }
#endif
        if (kernel_new::accepts(argv[1]))
            return kernel_new::main(argc, argv);
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        throw std::runtime_error("resident MPFR program accepts only new goal11 routine names");
#endif
#if defined(KERNEL_BENCH_MPFR_CUDA)
        // Small sizes too: the thresholds of the entry points do not apply.
        setenv("MPLAPACK_MPFR_CUDA_MIN_MNK", "0", 1);
        setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN", "0", 1);
        cupti_start();
#endif
#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_CUDA)
        int devices = 0;
        if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) {
            std::fprintf(stderr, "no CUDA device\n");
            return 77;
        }
        cudaFree(0);
#endif
        const spec s = parse(argv[1], argv[2]);
        const std::vector<long> sizes = list(argv[3]);
        return s.complex ? bench<COMPLEX>(s, sizes, argc, argv) : bench<REAL>(s, sizes, argc, argv);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
