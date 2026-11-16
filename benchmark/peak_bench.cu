// Arithmetic and memory peaks of the A100 for the multiple-precision types
// used by MPLAPACK's GPU kernels: the denominators of the roofline analysis
// (benchmark/roofline.py).
//
// For each type, the throughput (operations per second, whole GPU) of
//   madd: c = c + c*b   (the inner step of GEMM/SYRK/TRSM updates; the product changes every step)
//   add:  c = c + a
//   mul:  c = c * a     (a close to 1)
//   div:  c = a / c     (dd/td/qd: the library's division; MPFR fixed: cu_fdiv_rn)
// with the same arithmetic as the kernels (dd/td/qd: mplapack_gemm::{dd,td,qd}_ops,
// libQD3 QD_BF algorithms, compiled with -fmad=false; MPFR: cu_fp::cu_freal<PB>),
// each thread running ILP independent chains so that latency is hidden, and
// the device-memory bandwidth (read+write of a large array).
//
// Output (CSV): type,op,ops_per_second,note (the best over the launch shapes tried)
//
// build: nvcc -O3 -std=c++17 -fmad=false -arch=sm_80 -I../include -I../mpblas/optimized/common \
//        -I../external/i/QD/include -I../external/i/MPC_CUDA/include peak_bench.cu -o peak_bench

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "Rgemm_cuda_ops_common.h"
#include <mpc_cuda/cu_fround.cuh>

#define CHECK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { std::fprintf(stderr, "%s: %s\n", #x, cudaGetErrorString(e)); std::exit(1); } } while (0)

enum { MADD, ADD, MUL, DIV };
static const char *opname[] = {"madd", "add", "mul", "div"};

// ---- dd / td / qd ----
template <class Ops, int ILP> __global__ void bf_kernel(int op, long iters, typename Ops::F *out) {
    typedef typename Ops::F F;
    F c[ILP], a = Ops::from_d(1.0 + 1e-9 * threadIdx.x), b = Ops::from_d(1e-9);
    for (int i = 0; i < ILP; i++) c[i] = Ops::from_d(1.0 + 1e-7 * i);
    for (long it = 0; it < iters; it++) {
#pragma unroll
        for (int i = 0; i < ILP; i++) {
            if (op == MADD) c[i] = Ops::add(c[i], Ops::mul(c[i], b)); // the product changes every step: no hoisting
            else if (op == ADD) c[i] = Ops::add(c[i], a);
            else if (op == MUL) c[i] = Ops::mul(c[i], a);
            else c[i] = Ops::div(a, c[i]);
        }
    }
    F s = c[0];
    for (int i = 1; i < ILP; i++) s = Ops::add(s, c[i]);
    out[blockIdx.x * blockDim.x + threadIdx.x] = s;
}

// ---- MPFR fixed precision (cu_freal) ----
template <int PB, int ILP> __global__ void freal_kernel(int op, long iters, cu_fp::cu_freal<PB> *out) {
    typedef cu_fp::cu_freal<PB> F;
    F c[ILP], a(1.0 + 1e-9 * threadIdx.x), b(1e-9);
    for (int i = 0; i < ILP; i++) c[i] = F(1.0 + 1e-7 * i);
    for (long it = 0; it < iters; it++) {
        for (int i = 0; i < ILP; i++) {
            if (op == MADD) c[i] = c[i] + c[i] * b; // the product changes every step
            else if (op == ADD) c[i] = c[i] + a;
            else if (op == MUL) c[i] = c[i] * a;
            else { F q; cu_fp::cu_fdiv_rn<PB>(q, a, c[i]); c[i] = q; }
        }
    }
    F s = c[0];
    for (int i = 1; i < ILP; i++) s = s + c[i];
    out[blockIdx.x * blockDim.x + threadIdx.x] = s;
}

template <class K, class T> double run(K kernel, int op, long iters, int blocks, int threads) {
    T *out;
    CHECK(cudaMalloc(&out, sizeof(T) * (size_t)blocks * threads));
    cudaEvent_t e0, e1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    kernel<<<blocks, threads>>>(op, iters / 10 + 1, out); // warm-up
    CHECK(cudaDeviceSynchronize());
    double best = 1e30;
    for (int r = 0; r < 3; r++) {
        cudaEventRecord(e0);
        kernel<<<blocks, threads>>>(op, iters, out);
        cudaEventRecord(e1);
        CHECK(cudaEventSynchronize(e1));
        float ms;
        cudaEventElapsedTime(&ms, e0, e1);
        best = ms * 1e-3 < best ? ms * 1e-3 : best;
    }
    CHECK(cudaFree(out));
    return best;
}

// Best over thread-block sizes and occupancies: the peak, not one configuration.
template <class Ops, int ILP> double bf_one(int op, long it) {
    int sms;
    cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0);
    double best = 0;
    for (int threads : {64, 128, 256})
        for (int per_sm : {4, 8, 16}) {
            const int blocks = sms * per_sm;
            double t = run<decltype(&bf_kernel<Ops, ILP>), typename Ops::F>(bf_kernel<Ops, ILP>, op, it, blocks, threads);
            double r = (double)blocks * threads * ILP * it / t;
            best = r > best ? r : best;
        }
    return best;
}

template <class Ops> void bf(const char *name, long iters) {
    for (int op = MADD; op <= DIV; op++) {
        const long it = op == DIV ? iters / 8 : iters;
        double r = 0;
        for (double x : {bf_one<Ops, 1>(op, it), bf_one<Ops, 2>(op, it), bf_one<Ops, 4>(op, it)})
            r = x > r ? x : r;
        std::printf("%s,%s,%.6e,best of ILP 1/2/4, 64-256 threads\n", name, opname[op], r);
        std::fflush(stdout);
    }
}

template <int PB, int ILP> double freal_one(int op, long it) {
    int sms;
    cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0);
    double best = 0;
    for (int threads : {32, 64, 128})
        for (int per_sm : {4, 8, 16}) {
            const int blocks = sms * per_sm;
            double t = run<decltype(&freal_kernel<PB, ILP>), cu_fp::cu_freal<PB>>(freal_kernel<PB, ILP>, op, it, blocks, threads);
            double r = (double)blocks * threads * ILP * it / t;
            best = r > best ? r : best;
        }
    return best;
}

template <int PB> void freal(const char *name, long iters) {
    CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 16384));
    for (int op = MADD; op <= DIV; op++) {
        const long it = op == DIV ? iters / 8 : iters;
        double r = 0;
        for (double x : {freal_one<PB, 1>(op, it), freal_one<PB, 2>(op, it)})
            r = x > r ? x : r;
        std::printf("%s,%s,%.6e,best of ILP 1/2, 32-128 threads\n", name, opname[op], r);
        std::fflush(stdout);
    }
}

__global__ void copy_kernel(const double4 *__restrict__ in, double4 *__restrict__ out, size_t n) {
    for (size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x; i < n; i += (size_t)gridDim.x * blockDim.x)
        out[i] = in[i];
}

void bandwidth() {
    const size_t bytes = (size_t)4 << 30, n = bytes / sizeof(double4);
    double4 *a, *b;
    CHECK(cudaMalloc(&a, bytes));
    CHECK(cudaMalloc(&b, bytes));
    CHECK(cudaMemset(a, 0, bytes));
    cudaEvent_t e0, e1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    copy_kernel<<<108 * 16, 256>>>(a, b, n);
    CHECK(cudaDeviceSynchronize());
    double best = 1e30;
    for (int r = 0; r < 5; r++) {
        cudaEventRecord(e0);
        copy_kernel<<<108 * 16, 256>>>(a, b, n);
        cudaEventRecord(e1);
        CHECK(cudaEventSynchronize(e1));
        float ms;
        cudaEventElapsedTime(&ms, e0, e1);
        best = ms * 1e-3 < best ? ms * 1e-3 : best;
    }
    std::printf("memory,copy_bytes,%.6e,read+write of 2x4 GiB\n", 2.0 * bytes / best); // bytes read + written per second
    cudaFree(a);
    cudaFree(b);
}

int main(int argc, char **argv) {
    const long scale = argc > 1 ? std::atol(argv[1]) : 1;
    std::printf("type,op,ops_per_second,note\n");
    bandwidth();
    bf<mplapack_gemm::dd_ops>("dd", 20000 * scale);
    bf<mplapack_gemm::td_ops>("td", 6000 * scale);
    bf<mplapack_gemm::qd_ops>("qd", 3000 * scale);
    freal<256>("mpfr256", 2000 * scale);
    freal<512>("mpfr512", 800 * scale);
    freal<1024>("mpfr1024", 200 * scale);
    return 0;
}
