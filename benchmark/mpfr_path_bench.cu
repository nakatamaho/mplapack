// Path-by-path throughput of the MPFR fixed-precision arithmetic (cu_fp::cu_freal<PB>)
// on the GPU, and of the integer primitives it is made of: the inputs of the
// two-level MPFR roofline (benchmark/results/2026/<host>/roofline/mpfr/).
//
// Level 1 (hardware, implementation-independent): the throughput of
//   mad64   one 64x64 -> 128-bit limb product accumulated with carries
//           (mad.lo.cc / madc.hi / add.cc / addc: the inner step of a schoolbook product)
//   addc64  one 64-bit add with carry-in/out
//   shf64   one 64-bit funnel shift by a run-time amount ((x >> s) | (y << (64 - s)))
// so that a PB-bit multiplication needs at least N^2 mad64 (full product,
// N = PB/64 limbs; about N(N+1)/2 + N for a truncated product) and an addition
// N shf64 + N addc64.
//
// Level 2 (cu_freal as implemented): the throughput of cu_fmul and of cu_fadd
// on each of its paths, chosen by the operands:
//   mul            any (one path)
//   add_same_g*    same signs (cu_addmag), exponent gap g = 0, 1..63, 64..PB-1, >= PB
//   add_opp_g*     opposite signs (cu_submag), same gaps; g0 has no large cancellation
//   add_opp_cancel opposite signs, equal exponents, the top limbs equal: heavy cancellation
// Each thread keeps a fixed a and a b whose low limb changes every step (no
// hoisting); the gap and sign stay in their category.  Best over launch shapes.
//
// Output (CSV): pb,op,ops_per_second
// build: nvcc -O3 -std=c++17 -fmad=false -arch=sm_80 -U_GNU_SOURCE -D_XOPEN_SOURCE=700 \
//        -I../external/i/MPC_CUDA/include mpfr_path_bench.cu -o mpfr_path_bench

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <mpc_cuda/cu_freal.cuh>

#define CHECK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { std::fprintf(stderr, "%s: %s\n", #x, cudaGetErrorString(e)); std::exit(1); } } while (0)

typedef unsigned long long u64;

// ---- level 1 primitives ----
__global__ void k_mad64(long iters, u64 *out) {
    u64 a = 0x9e3779b97f4a7c15ULL ^ threadIdx.x, b = 0xbf58476d1ce4e5b9ULL, r0 = 1, r1 = 2, r2 = 3, r3 = 4, c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    for (long it = 0; it < iters; it++) {
        u64 lo, hi;
#define STEP(r, c)                                                                                         \
    asm volatile("mad.lo.cc.u64 %0,%2,%3,%4;\n\tmadc.hi.u64 %1,%2,%3,0;\n\tadd.cc.u64 %0,%0,%5;\n\taddc.u64 %1,%1,0;" \
                 : "=&l"(lo), "=&l"(hi) : "l"(a), "l"(b), "l"(c), "l"(r));                                   \
    r = lo; c = hi;
        STEP(r0, c0) STEP(r1, c1) STEP(r2, c2) STEP(r3, c3)
        a += 0x2545F4914F6CDD1DULL;
    }
    out[blockIdx.x * blockDim.x + threadIdx.x] = r0 ^ r1 ^ r2 ^ r3 ^ c0 ^ c1 ^ c2 ^ c3;
}
__global__ void k_addc64(long iters, u64 *out) {
    u64 x[8], y = 0x9e3779b97f4a7c15ULL ^ threadIdx.x;
    for (int i = 0; i < 8; i++) x[i] = i * 77 + threadIdx.x;
    for (long it = 0; it < iters; it++) {
        asm volatile("add.cc.u64 %0,%0,%8;\n\taddc.cc.u64 %1,%1,%8;\n\taddc.cc.u64 %2,%2,%8;\n\taddc.cc.u64 %3,%3,%8;\n\t"
                     "addc.cc.u64 %4,%4,%8;\n\taddc.cc.u64 %5,%5,%8;\n\taddc.cc.u64 %6,%6,%8;\n\taddc.u64 %7,%7,%8;"
                     : "+l"(x[0]), "+l"(x[1]), "+l"(x[2]), "+l"(x[3]), "+l"(x[4]), "+l"(x[5]), "+l"(x[6]), "+l"(x[7]) : "l"(y));
        y = y * 6364136223846793005ULL + 1;
    }
    u64 s = 0;
    for (int i = 0; i < 8; i++) s ^= x[i];
    out[blockIdx.x * blockDim.x + threadIdx.x] = s;
}
__global__ void k_shf64(long iters, int shift, u64 *out) {
    u64 x[8];
    for (int i = 0; i < 8; i++) x[i] = 0x9e3779b97f4a7c15ULL * (i + 1) ^ threadIdx.x;
    int s = shift + (threadIdx.x & 1);
    for (long it = 0; it < iters; it++) {
#pragma unroll
        for (int i = 0; i < 8; i++) x[i] = (x[i] >> s) | (x[(i + 1) & 7] << (64 - s));
        s = (s & 31) + 1;
    }
    u64 r = 0;
    for (int i = 0; i < 8; i++) r ^= x[i];
    out[blockIdx.x * blockDim.x + threadIdx.x] = r;
}

// ---- level 2: cu_freal paths ----
enum { MUL, ADD_SAME_G0, ADD_SAME_G1, ADD_SAME_G64, ADD_SAME_GBIG, ADD_OPP_G0, ADD_OPP_G1, ADD_OPP_G64, ADD_OPP_GBIG, ADD_OPP_CANCEL, NOPS };
static const char *opname[] = {"mul", "add_same_g0", "add_same_g1_63", "add_same_g64_pb", "add_same_gbig",
                               "add_opp_g0", "add_opp_g1_63", "add_opp_g64_pb", "add_opp_gbig", "add_opp_cancel"};

template <int PB> __global__ void k_freal(int op, long iters, unsigned long long *out) {
    typedef cu_fp::cu_freal<PB> F;
    const int N = F::N;
    F a, b;
    a.sign = 1;
    a.exp = 0;
    for (int i = 0; i < N; i++) a.m[i] = 0x9e3779b97f4a7c15ULL * (i + 3) ^ (threadIdx.x * 0x100000001ULL);
    a.m[N - 1] |= 1ULL << 63;
    if (F::SB) a.m[0] &= ~0xffffffffULL;
    b = a;
    long gap = 0;
    switch (op) {
    case ADD_SAME_G1: case ADD_OPP_G1: gap = 17; break;
    case ADD_SAME_G64: case ADD_OPP_G64: gap = PB / 2 + 5; break;
    case ADD_SAME_GBIG: case ADD_OPP_GBIG: gap = PB + 40; break;
    default: gap = 0;
    }
    b.exp = a.exp - gap;
    if (op >= ADD_OPP_G0) b.sign = -1;
    if (op == ADD_OPP_G0) b.m[N - 1] = (b.m[N - 1] >> 1) | (1ULL << 63);      // |b| < |a|, top limb differs: little cancellation
    if (op == ADD_OPP_CANCEL) b.m[0] ^= 0x10000ULL;                          // only the low limb differs: heavy cancellation
    u64 acc = 0;
    for (long it = 0; it < iters; it++) {
        F x;
        u64 low = (u64)(it & 0xff) << (F::SB + 1);
        b.m[0] = (b.m[0] & ~(0xffULL << (F::SB + 1))) | low; // the low limb changes; gap/sign category unchanged
        if (op == MUL)
            x = a * b;
        else
            x = a + b;
        acc ^= x.m[N - 1] ^ (u64)x.exp;
    }
    out[blockIdx.x * blockDim.x + threadIdx.x] = acc;
}

template <class K, class... A> double best_rate(K kernel, long work_per_thread, A... args) {
    int sms;
    cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0);
    u64 *out;
    CHECK(cudaMalloc(&out, sizeof(u64) * sms * 16 * 256));
    cudaEvent_t e0, e1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    double best = 0;
    for (int threads : {64, 128, 256})
        for (int per_sm : {2, 4, 8, 16}) {
            const int blocks = sms * per_sm;
            kernel<<<blocks, threads>>>(args..., out); // warm-up
            CHECK(cudaDeviceSynchronize());
            for (int r = 0; r < 2; r++) {
                cudaEventRecord(e0);
                kernel<<<blocks, threads>>>(args..., out);
                cudaEventRecord(e1);
                CHECK(cudaEventSynchronize(e1));
                float ms;
                cudaEventElapsedTime(&ms, e0, e1);
                double rate = (double)blocks * threads * work_per_thread / (ms * 1e-3);
                best = rate > best ? rate : best;
            }
        }
    CHECK(cudaFree(out));
    return best;
}

template <int PB> void freal(long iters) {
    CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 16384));
    for (int op = 0; op < NOPS; op++) {
        long it = op == MUL ? iters : iters * 2;
        double r = best_rate(k_freal<PB>, it, op, it);
        std::printf("%d,%s,%.6e\n", PB, opname[op], r);
        std::fflush(stdout);
    }
}

int main(int argc, char **argv) {
    const long scale = argc > 1 ? std::atol(argv[1]) : 1;
    std::printf("pb,op,ops_per_second\n");
    std::printf("0,mad64,%.6e\n", best_rate(k_mad64, 4 * 20000 * scale, 20000L * scale));
    std::printf("0,addc64,%.6e\n", best_rate(k_addc64, 8 * 20000 * scale, 20000L * scale));
    std::printf("0,shf64,%.6e\n", best_rate(k_shf64, 8 * 20000 * scale, 20000L * scale, 7));
    std::fflush(stdout);
    freal<256>(400 * scale);
    freal<512>(200 * scale);
    freal<1024>(60 * scale);
    return 0;
}
