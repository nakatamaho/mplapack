// Internal device building blocks; compile this TU with C++ (emulation)
// or nvcc. Oracle panels are the unmodified FABLE reference bodies.
#if defined(__CUDACC__)
// nvcc already provides the two-argument long min/max overloads. These
// reference panel bodies need no Fortran variadic overloads.
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#endif
#include <mpblas.h>
#include <mplapack.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <limits>
#include "../Rlapack_cuda_common.cuh"
#define Rgetrf2 panel_reference_getrf2
#include "../../../../mplapack/reference/Rgetrf2.cpp"
#undef Rgetrf2
#define Rpotrf2 panel_reference_potrf2
#include "../../../../mplapack/reference/Rpotrf2.cpp"
#undef Rpotrf2

namespace {
using namespace mplapack_lapack_gpu;
#if defined(MPLAPACK_BUILD_WITH_DD)
using Ops = mplapack_gemm::dd_ops;
#elif defined(MPLAPACK_BUILD_WITH_TD)
using Ops = mplapack_gemm::td_ops;
#elif defined(MPLAPACK_BUILD_WITH_QD)
using Ops = mplapack_gemm::qd_ops;
#else
#error dd/td/qd only
#endif
using R = dreal<Ops>;
using C = dcomplex<R>;
long checks = 0;
void require(bool ok, const char *msg) {
    ++checks;
    if (!ok) {
        std::cerr << "FAIL: " << msg << " check=" << checks << "\n";
        std::exit(1);
    }
}
void success(int e) {
    if (e)
        std::cerr << "runtime error " << e << "\n";
    require(e == 0, "runtime");
}
R pack(const REAL &a) {
    R r;
    for (unsigned q = 0; q < sizeof(r.v.x) / sizeof(double); ++q)
        r.v.x[q] = a.x[q];
    return r;
}
REAL unpack(const R &a) {
    REAL r;
    for (unsigned q = 0; q < sizeof(a.v.x) / sizeof(double); ++q)
        r.x[q] = a.v.x[q];
    return r;
}
void equal(const std::vector<R> &a, const std::vector<REAL> &b) {
    require(a.size() == b.size(), "size");
    for (size_t i = 0; i < a.size(); ++i) {
        R r = pack(b[i]);
        if (std::memcmp(&a[i], &r, sizeof(R))) {
            std::cerr << "element " << i << "\n";
            require(false, "bit identity");
        }
    }
    ++checks;
}
std::vector<R> packed(const std::vector<REAL> &a) {
    std::vector<R> p;
    for (const auto &v : a)
        p.push_back(pack(v));
    return p;
}
REAL value(long i) {
    return REAL((i * 37 % 101) - 50) / REAL(13.0) + REAL((i * 17 % 19) - 9) * REAL(0x1p-80) + REAL((i * 11 % 13) - 6) * REAL(0x1p-160) +
           REAL((i * 3 % 7) - 3) * REAL(0x1p-240);
}
struct roots {
    R *a;
    MPLAPACK_DEV_HD void operator()(long i) const { a[i] = R(square_root<Ops>(a[i].v)); }
};
void test_allocation() {
    buffer<R> b;
    require(b.allocate(std::numeric_limits<size_t>::max()) != 0, "allocation byte overflow");
    matrix<R> a;
    require(a.allocate(std::numeric_limits<long>::max(), std::numeric_limits<long>::max()) != 0, "matrix size overflow");
    panel<Ops> bad(-1, 1);
    require(bad.error != 0, "negative workspace");
}
void test_sqrt() {
    std::vector<REAL> a = {REAL(0.0), REAL(-0.0), REAL(1.0), REAL(2.0), REAL(0x1p+970), REAL(0x1p+1020), REAL(0x1p-970)};
    a.push_back(REAL(std::numeric_limits<double>::denorm_min()));
    a.push_back(REAL(std::numeric_limits<double>::max()));
    a.push_back(REAL(std::numeric_limits<double>::infinity()));
    a.push_back(REAL(std::numeric_limits<double>::quiet_NaN()));
    a.push_back(REAL(-1.0));
    for (long i = 1; i <= 256; ++i)
        a.push_back(abs(value(i)) + REAL(0x1p-100));
    for (int e = -1000; e <= 1000; e += 37)
        a.push_back((REAL(1.0) + REAL(0x1p-54)) * REAL(std::ldexp(1.0, e)));
    auto p = packed(a);
    matrix<R> d;
    success(d.allocate(p.size(), 1));
    success(d.copy(p.data(), p.size(), false));
    success(elements(roots{d.data()}, p.size()));
    success(d.copy(p.data(), p.size(), true));
    for (auto &v : a)
        v = sqrt(v);
    equal(p, a);
}
void test_pivots() {
    for (long n : {0L, 1L, 2L, 3L, 255L, 256L, 257L, 513L, 1025L}) {
        std::vector<REAL> a(std::max(1L, n));
        for (long i = 0; i < n; ++i)
            a[i] = value(i);
        if (n > 2) {
            a[n - 1] = REAL(-100.0);
            a[n / 2] = REAL(100.0);
        }
        auto p = packed(a);
        buffer<R> d;
        buffer<long> idx;
        success(d.allocate(p.size()));
        success(idx.allocate(1));
        success(transfer(d.data(), p.data(), p.size(), false));
        success(iamax(d.data(), n, idx.data()));
        long got = -1;
        success(transfer(&got, idx.data(), 1, true));
        require(got == iRamax(n, a.data(), 1), "iRamax tie/length");
        std::vector<COMPLEX> ca(a.size());
        std::vector<C> cp(a.size());
        for (long i = 0; i < n; ++i) {
            ca[i] = COMPLEX(a[i], a[i] / REAL(2.0));
            cp[i] = C(pack(ca[i].real()), pack(ca[i].imag()));
        }
        buffer<C> cd;
        success(cd.allocate(cp.size()));
        success(transfer(cd.data(), cp.data(), cp.size(), false));
        success(iamax(cd.data(), n, idx.data()));
        success(transfer(&got, idx.data(), 1, true));
        require(got == iCamax(n, ca.data(), 1), "iCamax");
        if (n > 2) {
            for (long pos : {0L, n / 2}) {
                a[pos] = REAL(std::numeric_limits<double>::quiet_NaN());
                p = packed(a);
                success(transfer(d.data(), p.data(), p.size(), false));
                success(iamax(d.data(), n, idx.data()));
                success(transfer(&got, idx.data(), 1, true));
                require(got == iRamax(n, a.data(), 1), "iRamax NaN");
                a[pos] = value(pos);
            }
        }
    }
}
void test_scaling() {
    const REAL sfmin = Rlamch("S");
    for (REAL pivot : {REAL(2.0), REAL(-3.0), value(5), sfmin, sfmin / REAL(2.0)}) {
        std::vector<REAL> a = {pivot, value(7) * pivot / REAL(100.0), -pivot / REAL(7.0), REAL(0.0)};
        auto p = packed(a);
        matrix<R> d;
        success(d.allocate(4, 1));
        success(d.copy(p.data(), 4, false));
        buffer<R> reciprocal;
        success(reciprocal.allocate(1));
        success(scale<Ops>(d.data(), 4, pack(sfmin), reciprocal.data()));
        success(d.copy(p.data(), 4, true));
        if (abs(pivot) >= sfmin)
            Rscal(3, REAL(1.0) / pivot, a.data() + 1, 1);
        else
            for (int i = 1; i < 4; ++i)
                a[i] = a[i] / pivot;
        equal(p, a);
    }
}
void test_swaps() {
    const long m = 9, n = 37, ld = 13;
    for (long inc : {1L, -1L, 2L, -2L, 0L}) {
        std::vector<REAL> a(ld * n);
        for (long i = 0; i < ld * n; ++i)
            a[i] = value(i);
        auto p = packed(a);
        long piv[18];
        INTEGER refp[18];
        for (int i = 0; i < 18; ++i)
            piv[i] = refp[i] = 1 + (i * 5 % 9);
        matrix<R> d;
        buffer<long> dp;
        success(d.allocate(m, n));
        success(dp.allocate(18));
        success(d.copy(p.data(), ld, false));
        success(transfer(dp.data(), piv, 18, false));
        success(laswp(n, d.data(), d.ld, 2, 7, dp.data(), inc));
        success(d.copy(p.data(), ld, true));
        Rlaswp(n, a.data(), ld, 2, 7, refp, inc);
        equal(p, a);
    }
}
void test_lu(long m, long n, int kind) {
    long ld = m + 3;
    std::vector<REAL> a(ld * n, REAL(-777.0));
    for (long j = 0; j < n; ++j)
        for (long i = 0; i < m; ++i)
            a[i + j * ld] = value(i + j * m + 1);
    if (kind == 1)
        for (long i = 0; i < m; ++i)
            a[i + (n / 2) * ld] = REAL(0.0);
    if (kind == 2)
        for (long j = 0; j < n; ++j)
            for (long i = 0; i < m; ++i)
                a[i + j * ld] = REAL((i + j) % 2 ? -1.0 : 1.0);
    if (kind == 4)
        for (long j = 0; j < n; ++j)
            for (long i = 0; i < m; ++i)
                a[i + j * ld] = i == j + 1 ? REAL(m + 1) : value(i + j + 1) / REAL(100.0);
    if (kind == 3)
        for (long j = 0; j < n; ++j)
            for (long i = 0; i < m; ++i)
                a[i + j * ld] *= Rlamch("S") / REAL(4.0);
    const auto original = a;
    auto p = packed(a);
    matrix<R> d;
    buffer<long> piv;
    success(d.allocate(m, n));
    success(d.copy(p.data(), ld, false));
    success(piv.allocate(std::min(m, n)));
    panel<Ops> ctx(m, n);
    success(ctx.error);
    long info = ctx.getrf2(m, n, d.data(), d.ld, piv.data(), pack(Rlamch("S")));
    success(ctx.error);
    success(d.copy(p.data(), ld, true));
    std::vector<long> ip(std::min(m, n));
    success(transfer(ip.data(), piv.data(), ip.size(), true));
    INTEGER ri = 0;
    std::vector<INTEGER> rp(ip.size());
    panel_reference_getrf2(m, n, a.data(), ld, rp.data(), ri);
    require(info == ri, "LU info");
    for (size_t i = 0; i < ip.size(); ++i)
        require(ip[i] == rp[i], "LU pivots");
    equal(p, a);
    if (kind == 0 && m <= 33 && n <= 33) {
        auto pa = original;
        for (long k = 0; k < std::min(m, n); ++k)
            for (long j = 0; j < n; ++j)
                std::swap(pa[k + j * ld], pa[ip[k] - 1 + j * ld]);
        REAL norm = 0.0, err = 0.0;
        for (long j = 0; j < n; ++j) {
            REAL asum = 0.0, rsum = 0.0;
            for (long i = 0; i < m; ++i) {
                REAL sum = 0.0;
                for (long k = 0; k < std::min(m, n); ++k)
                    if (i >= k && k <= j)
                        sum += (i == k ? REAL(1.0) : unpack(p[i + k * ld])) * unpack(p[k + j * ld]);
                asum += abs(original[i + j * ld]);
                rsum += abs(pa[i + j * ld] - sum);
            }
            if (asum > norm)
                norm = asum;
            if (rsum > err)
                err = rsum;
        }
        require(norm == REAL(0.0) || err / (REAL(n) * norm * Rlamch("E")) <= REAL(30.0), "LU backward error <=30");
    }
}
void test_chol(long n, bool upper, int bad) {
    long ld = n + 3;
    std::vector<REAL> a(ld * n, REAL(-777.0));
    // Strictly diagonally dominant SPD, with nontrivial low words.
    for (long j = 0; j < n; ++j)
        for (long i = 0; i < n; ++i)
            if (upper ? i <= j : i >= j)
                a[i + j * ld] = i == j ? REAL(n * 10 + 3) : value(i + j + 1) / REAL(100.0);
    if (bad)
        a[(bad - 1) * (ld + 1)] = REAL(-1.0);
    const auto original = a;
    auto p = packed(a);
    matrix<R> d;
    success(d.allocate(n, n));
    success(d.copy(p.data(), ld, false));
    panel<Ops> ctx(n, n);
    success(ctx.error);
    long info = ctx.potrf2(upper, n, d.data(), d.ld);
    success(ctx.error);
    success(d.copy(p.data(), ld, true));
    INTEGER ri = 0;
    panel_reference_potrf2(upper ? "U" : "L", n, a.data(), ld, ri);
    require(info == ri, "Cholesky info");
    equal(p, a);
    if (!bad && n > 0 && n <= 33) {
        REAL norm = 0.0, err = 0.0;
        for (long j = 0; j < n; ++j) {
            REAL asum = 0.0, rsum = 0.0;
            for (long i = 0; i < n; ++i) {
                const long r = upper ? std::min(i, j) : std::max(i, j), c = upper ? std::max(i, j) : std::min(i, j);
                const REAL orig = original[r + c * ld];
                REAL sum = 0.0;
                for (long k = 0; k <= std::min(i, j); ++k)
                    sum += unpack(p[upper ? k + i * ld : i + k * ld]) * unpack(p[upper ? k + j * ld : j + k * ld]);
                asum += abs(orig);
                rsum += abs(orig - sum);
            }
            if (asum > norm)
                norm = asum;
            if (rsum > err)
                err = rsum;
        }
        require(err / (REAL(n) * norm * Rlamch("E")) <= REAL(30.0), "Cholesky backward error <=30");
    }
}
void test_herk() {
    using J = job_Cherk<Ops>;
    for (bool upper : {false, true})
        for (bool trans : {false, true}) {
            const long n = 7, k = 5, ld = 11;
            std::vector<COMPLEX> a(ld * std::max(n, k)), c(ld * n);
            std::vector<C> pa(a.size()), pc(c.size());
            for (size_t i = 0; i < a.size(); ++i) {
                a[i] = COMPLEX(value(i), value(i + 7));
                pa[i] = C(pack(a[i].real()), pack(a[i].imag()));
            }
            for (size_t i = 0; i < c.size(); ++i) {
                c[i] = COMPLEX(value(i + 3), value(i + 9));
                pc[i] = C(pack(c[i].real()), pack(c[i].imag()));
            }
            buffer<C> da, dc;
            buffer<J> scratch;
            success(da.allocate(pa.size()));
            success(dc.allocate(pc.size()));
            success(scratch.allocate(1));
            success(transfer(da.data(), pa.data(), pa.size(), false));
            success(transfer(dc.data(), pc.data(), pc.size(), false));
            J job{};
            job.uplo[0] = upper ? 'U' : 'L';
            job.tr[0] = trans ? 'C' : 'N';
            job.n = n;
            job.k = k;
            job.lda = ld;
            job.ldc = ld;
            job.alpha = R(-1.0);
            job.beta = R(1.0);
            job.A = da.data();
            job.C = dc.data();
            success(job_device(job, job.units(), scratch.data()));
            success(transfer(pc.data(), dc.data(), pc.size(), true));
            Cherk(upper ? "U" : "L", trans ? "C" : "N", n, k, REAL(-1.0), a.data(), ld, REAL(1.0), c.data(), ld);
            for (size_t i = 0; i < c.size(); ++i) {
                C expected(pack(c[i].real()), pack(c[i].imag()));
                require(!std::memcmp(&pc[i], &expected, sizeof(C)), "Cherk device pointers");
            }
        }
}
} // namespace
int main() {
#if defined(__CUDACC__)
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0)
        return 77;
#endif
    test_allocation();
    std::cout << "sqrt\n";
    test_sqrt();
    std::cout << "pivots\n";
    test_pivots();
    std::cout << "scaling\n";
    test_scaling();
    std::cout << "swaps\n";
    test_swaps();
    std::cout << "herk\n";
    test_herk();
    std::cout << "panels\n";
    for (long n : {1L, 2L, 3L, 15L, 16L, 17L, 31L, 32L, 33L, 63L, 64L, 65L, 67L, 131L}) {
        for (int kind = 0; kind < 4; ++kind) {
            test_lu(n, n, kind);
            test_lu(n + 3, n, kind);
            test_lu(n, n + 3, kind);
        }
        for (bool upper : {false, true}) {
            test_chol(n, upper, 0);
            test_chol(n, upper, 1);
            test_chol(n, upper, n);
            test_chol(n, upper, (n + 1) / 2);
        }
    }
    test_lu(0, 3, 0);
    test_lu(3, 0, 0);
    test_chol(0, false, 0);
    test_lu(83, 67, 4);
    test_lu(599, 33, 0);
    test_lu(33, 599, 0);
    test_chol(259, false, 0);
    test_chol(259, true, 0);
    test_chol(599, false, 0);
    test_chol(599, true, 0);
    std::cout << "PASS: device LAPACK blocks, " << checks << " checks\n";
}
