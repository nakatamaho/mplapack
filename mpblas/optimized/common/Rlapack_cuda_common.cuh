/* Internal dd/td/qd device-resident LAPACK building blocks.
 * No public entry points. All arithmetic uses the existing BF device BLAS.
 * Calls are synchronous and return CUDA error codes (0 in host emulation).
 * A failed internal operation may change device buffers; callers must not
 * publish/download results until the entire operation succeeds.
 */
#ifndef MPLAPACK_LAPACK_CUDA_COMMON_CUH
#define MPLAPACK_LAPACK_CUDA_COMMON_CUH
#include <limits>
#include <algorithm>
#include <new>
#include "Rgemm_cuda_common.cuh"
#include "Rblas_cuda_common.cuh"
#include "Rblas_dev_types.h"
#include "Cblas_cuda_common.cuh"
#include "Rlapack_sqrt_common.h"

namespace mplapack_lapack_gpu {
using namespace mplapack_blas_dev;

// Owning, noncopyable storage. Overflow and available memory are checked
// before allocation. The matrix plus caller-owned scratch is the full budget.
template <class E> class buffer {
    E *p_ = nullptr;

  public:
    buffer() = default;
    buffer(const buffer &) = delete;
    buffer &operator=(const buffer &) = delete;
    ~buffer() {
#if defined(__CUDACC__)
        cudaFree(p_);
#else
        delete[] p_;
#endif
    }
    E *data() const { return p_; }
    int allocate(size_t count) {
        if (p_ || count > std::numeric_limits<size_t>::max() / sizeof(E))
            return 1;
        count = std::max(size_t(1), count);
#if defined(__CUDACC__)
        size_t free = 0, total = 0;
        cudaError_t e = cudaMemGetInfo(&free, &total);
        if (e != cudaSuccess)
            return (int)e;
        if (count * sizeof(E) > free)
            return (int)cudaErrorMemoryAllocation;
        return (int)cudaMalloc(&p_, count * sizeof(E));
#else
        p_ = new (std::nothrow) E[count];
        return p_ ? 0 : 1;
#endif
    }
};
template <class E> static int transfer(E *dst, const E *src, size_t n, bool down) {
    if (!n)
        return 0;
#if defined(__CUDACC__)
    return (int)cudaMemcpy(dst, src, n * sizeof(E), down ? cudaMemcpyDeviceToHost : cudaMemcpyHostToDevice);
#else
    (void)down;
    std::copy(src, src + n, dst);
    return 0;
#endif
}
// Packed column-major device storage, ld=rows. Padding in the host matrix
// is neither read nor written. Also usable with device real/complex elements.
template <class E> class matrix {
    buffer<E> storage;

  public:
    long rows = 0, cols = 0, ld = 1;
    E *data() const { return storage.data(); }
    int allocate(long m, long n) {
        if (m < 0 || n < 0)
            return 1;
        const size_t stride = (size_t)std::max(1L, m);
        if ((size_t)n > std::numeric_limits<size_t>::max() / stride)
            return 1;
        rows = m;
        cols = n;
        ld = (long)stride;
        return storage.allocate(stride * (size_t)n);
    }
    int copy(E *host, long lda, bool down) {
        if (lda < std::max(1L, rows))
            return 1;
        if (!rows || !cols)
            return 0;
#if defined(__CUDACC__)
        return (int)cudaMemcpy2D(down ? host : data(), (down ? lda : ld) * sizeof(E), down ? data() : host, (down ? ld : lda) * sizeof(E), rows * sizeof(E),
                                 cols, down ? cudaMemcpyDeviceToHost : cudaMemcpyHostToDevice);
#else
        for (long j = 0; j < cols; ++j)
            if (down)
                std::copy(data() + j * ld, data() + j * ld + rows, host + j * lda);
            else
                std::copy(host + j * lda, host + j * lda + rows, data() + j * ld);
        return 0;
#endif
    }
};

#if defined(__CUDACC__)
template <class Job> __global__ void k_elements(Job job, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        job(i);
}
#endif
template <class Job> static int elements(Job job, long n) {
    if (n <= 0)
        return 0;
#if defined(__CUDACC__)
    k_elements<<<(unsigned)((n + 63) / 64), 64>>>(job, n);
    cudaError_t e = cudaGetLastError();
    return (int)(e == cudaSuccess ? cudaDeviceSynchronize() : e);
#else
    for (long i = 0; i < n; ++i)
        job(i);
    return 0;
#endif
}

template <class Ops> MPLAPACK_DEV_HD dreal<Ops> magnitude(const dreal<Ops> &a) { return abs(a); }
template <class R> MPLAPACK_DEV_HD R magnitude(const dcomplex<R> &a) { return RCabs1(a); }

// One block scans an arbitrarily long column, then reduces 256 candidates.
// Index zero is sticky if NaN, matching the serial iRamax/iCamax test;
// other NaNs are ignored. Equal magnitudes choose the first row.
template <class E> MPLAPACK_DEV_HD long pivot_part(const E *a, long n, long start, long step) {
    long best = -1;
    for (long i = start; i < n; i += step) {
        const auto v = magnitude(a[i]);
        if (i == 0 || v > decltype(v)(-1.0)) {
            if (best < 0 || magnitude(a[i]) > magnitude(a[best]))
                best = i;
        }
    }
    return best;
}
template <class E> MPLAPACK_DEV_HD long pivot_merge(const E *a, long x, long y) {
    if (x < 0)
        return y;
    if (y < 0)
        return x;
    if (x > y) {
        long t = x;
        x = y;
        y = t;
    }
    return magnitude(a[y]) > magnitude(a[x]) ? y : x;
}
#if defined(__CUDACC__)
template <class E> __global__ void k_iamax(const E *a, long n, long *out) {
    __shared__ long idx[256];
    int t = threadIdx.x;
    idx[t] = pivot_part(a, n, t, 256);
    __syncthreads();
    for (int s = 128; s; s /= 2) {
        if (t < s)
            idx[t] = pivot_merge(a, idx[t], idx[t + s]);
        __syncthreads();
    }
    if (t == 0)
        *out = n > 0 ? idx[0] + 1 : 0;
}
#endif
template <class E> static int iamax(const E *a, long n, long *device_index) {
#if defined(__CUDACC__)
    k_iamax<<<1, 256>>>(a, n, device_index);
    cudaError_t e = cudaGetLastError();
    return (int)(e == cudaSuccess ? cudaDeviceSynchronize() : e);
#else
    long idx[256];
    for (int i = 0; i < 256; ++i)
        idx[i] = pivot_part(a, n, i, 256);
    for (int s = 128; s; s /= 2)
        for (int i = 0; i < s; ++i)
            idx[i] = pivot_merge(a, idx[i], idx[i + s]);
    *device_index = n > 0 ? idx[0] + 1 : 0;
    return 0;
#endif
}

template <class E> struct swaps {
    E *a;
    long ld, k1, k2, inc;
    const long *piv;
    MPLAPACK_DEV_HD void operator()(long j) const {
        long ix = inc > 0 ? k1 - 1 : k1 - 1 + (k1 - k2) * inc;
        for (long k = inc > 0 ? k1 : k2; inc > 0 ? k <= k2 : k >= k1; k += inc > 0 ? 1 : -1, ix += inc) {
            long p = piv[ix] - 1;
            E tmp = a[k - 1 + j * ld];
            a[k - 1 + j * ld] = a[p + j * ld];
            a[p + j * ld] = tmp;
        }
    }
};
// LAPACK 1-based pivots, including reverse/strided pivot traversal.
template <class E> static int laswp(long n, E *a, long ld, long k1, long k2, const long *piv, long inc) {
    return inc == 0 ? 0 : elements(swaps<E>{a, ld, k1, k2, inc, piv}, n);
}
template <class Ops> struct prepare_scale {
    using R = dreal<Ops>;
    R *a, *reciprocal;
    R sfmin;
    MPLAPACK_DEV_HD void operator()(long) const {
        if (a[0] != R(0.0) && abs(a[0]) >= sfmin)
            *reciprocal = R(1.0) / a[0];
    }
};
template <class Ops> struct scale_column {
    using R = dreal<Ops>;
    R *a;
    R sfmin;
    const R *reciprocal;
    MPLAPACK_DEV_HD void operator()(long i) const {
        // Swapping the pivot is a separate launch, so a[0] is stable here.
        if (a[0] != R(0.0)) {
            if (abs(a[0]) >= sfmin)
                a[i + 1] = *reciprocal * a[i + 1]; // Rscal: reciprocal is the left operand
            else
                a[i + 1] = a[i + 1] / a[0];
        }
    }
};
template <class Ops> static int scale(dreal<Ops> *a, long m, dreal<Ops> sfmin, dreal<Ops> *reciprocal) {
    if (m <= 1) return 0;
    int error = elements(prepare_scale<Ops>{a, reciprocal, sfmin}, 1);
    return error ? error : elements(scale_column<Ops>{a, sfmin, reciprocal}, m - 1);
}
template <class Ops> struct base_lu {
    using R = dreal<Ops>;
    R *a;
    long *piv, *info;
    bool single_row;
    MPLAPACK_DEV_HD void operator()(long) const {
        if (single_row)
            *piv = 1;
        const long p = *piv - 1;
        *info = a[p] == R(0.0) ? 1 : 0;
        if (!*info && p != 0) {
            R t = a[0];
            a[0] = a[p];
            a[p] = t;
        }
    }
};
template <class Ops> struct base_chol {
    using R = dreal<Ops>;
    R *a;
    long *info;
    MPLAPACK_DEV_HD void operator()(long) const {
        *info = (a[0] <= R(0.0) || Ops::is_nan(a[0].v)) ? 1 : 0;
        if (!*info)
            a[0] = R(square_root<Ops>(a[0].v));
    }
};
struct adjust_pivots {
    long *p;
    long offset;
    MPLAPACK_DEV_HD void operator()(long i) const { p[i] += offset; }
};

// Scratch is explicit and reusable. Device scalar feedback is limited to a
// pivot/status word; no panel entries are downloaded or computed on the CPU.
template <class Ops> class panel {
  public:
    using F = typename Ops::F;
    using R = dreal<Ops>;
    buffer<F> work, constants;
    buffer<long> status;
    int error = 0;
    size_t capacity = 0;
    panel(long m, long n) {
        if (m < 0 || n < 0 || (n && (size_t)m > std::numeric_limits<size_t>::max() / (size_t)n)) {
            error = 1;
            return;
        }
        capacity = (size_t)m * (size_t)n;
        error = work.allocate(capacity);
        if (!error)
            error = constants.allocate(3);
        if (!error)
            error = status.allocate(1);
        F c[2] = {Ops::from_d(1.0), Ops::neg(Ops::from_d(1.0))};
        if (!error)
            error = transfer(constants.data(), c, 2, false);
    }
    int trsm(bool left, bool upper, bool trans, bool unit, long m, long n, R *a, long lda, R *b, long ldb) {
        mplapack_blas::gpu_lane_args<Ops> g{};
        if (error)
            return error;
        g.s.op = mplapack_blas::L3_TRSM;
        g.s.left = left;
        g.s.upper = upper;
        g.s.trans = trans;
        g.s.nounit = !unit;
        g.s.m = m;
        g.s.n = n;
        g.s.alpha_one = 1;
        g.a = reinterpret_cast<F *>(a);
        g.c = reinterpret_cast<F *>(b);
        g.al = constants.data();
        g.lda = lda;
        g.ldc = ldb;
        return error = mplapack_blas::blas_device<Ops>(g);
    }
    int syrk(bool upper, bool trans, long n, long k, R *a, long lda, R *c, long ldc) {
        mplapack_blas::gpu_lane_args<Ops> g{};
        if (error)
            return error;
        g.s.op = mplapack_blas::L3_SYRK;
        g.s.upper = upper;
        g.s.trans = trans;
        g.s.n = n;
        g.s.k = k;
        g.s.beta_one = 1;
        g.a = reinterpret_cast<F *>(a);
        g.c = reinterpret_cast<F *>(c);
        g.al = constants.data() + 1;
        g.be = constants.data();
        g.lda = lda;
        g.ldc = ldc;
        return error = mplapack_blas::blas_device<Ops>(g);
    }
    int gemm(long m, long n, long k, R *a, long lda, R *b, long ldb, R *c, long ldc) {
        if (error)
            return error;
        if (m < 0 || n < 0 || k < 0 || (k && (size_t)n > capacity / (size_t)k))
            return error = 1;
        mplapack_gemm::shape s{};
        s.m = m;
        s.n = n;
        s.k = k;
        s.lda = lda;
        s.ldb = ldb;
        s.ldc = ldc;
        s.nota = true;
        s.notb = true;
        mplapack_gemm::gpu_scalars<F> sc{Ops::neg(Ops::from_d(1.0)), Ops::from_d(1.0), false, true};
        return error = mplapack_gemm::gemm_device<Ops>(s, sc, reinterpret_cast<F *>(a), reinterpret_cast<F *>(b), reinterpret_cast<F *>(c), work.data());
    }
    long getrf2(long m, long n, R *a, long ld, long *piv, R sfmin) {
        if (error || m == 0 || n == 0)
            return 0;
        if (m == 1 || n == 1) {
            if (m != 1)
                error = iamax(a, m, piv);
            if (!error)
                error = elements(base_lu<Ops>{a, piv, status.data(), m == 1}, 1);
            long info = 0;
            if (!error)
                error = transfer(&info, status.data(), 1, true);
            if (!error && !info && m > 1)
                error = scale<Ops>(a, m, sfmin, reinterpret_cast<R *>(constants.data() + 2));
            return info;
        }
        long n1 = std::min(m, n) / 2, n2 = n - n1;
        long info = getrf2(m, n1, a, ld, piv, sfmin);
        if (error)
            return info;
        error = laswp(n2, a + n1 * ld, ld, 1, n1, piv, 1);
        if (!error)
            trsm(true, false, false, true, n1, n2, a, ld, a + n1 * ld, ld);
        if (!error)
            gemm(m - n1, n2, n1, a + n1, ld, a + n1 * ld, ld, a + n1 + n1 * ld, ld);
        long i = getrf2(m - n1, n2, a + n1 + n1 * ld, ld, piv + n1, sfmin);
        if (!info && i)
            info = i + n1;
        if (!error)
            error = elements(adjust_pivots{piv + n1, n1}, std::min(m, n) - n1);
        if (!error)
            error = laswp(n1, a, ld, n1 + 1, std::min(m, n), piv, 1);
        return info;
    }
    long potrf2(bool upper, long n, R *a, long ld) {
        if (error || !n)
            return 0;
        if (n == 1) {
            error = elements(base_chol<Ops>{a, status.data()}, 1);
            long info = 0;
            if (!error)
                error = transfer(&info, status.data(), 1, true);
            return info;
        }
        long n1 = n / 2, n2 = n - n1;
        long info = potrf2(upper, n1, a, ld);
        if (error || info)
            return info;
        R *off = upper ? a + n1 * ld : a + n1;
        trsm(upper, upper, true, false, upper ? n1 : n2, upper ? n2 : n1, a, ld, off, ld);
        if (!error)
            syrk(upper, upper, n2, n1, off, ld, a + n1 + n1 * ld, ld);
        info = potrf2(upper, n2, a + n1 + n1 * ld, ld);
        return info ? info + n1 : 0;
    }
};
} // namespace mplapack_lapack_gpu
#endif
