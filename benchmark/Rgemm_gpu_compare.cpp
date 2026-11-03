// Times Rgemm/Cgemm("N","N",n,n,n,1,A,n,B,n,0,C,n), entries uniform in [-1,1], best of R runs.
// usage: gemm_bench R|C n [reps] [dumpfile]   -> best seconds, GFLOPS (Rgemm 2n^3/t, Cgemm 8n^3/t)
//        gemm_bench cmp R|C n ref other       -> max|ref-other| / max|ref|
// MPFR: precision from MPFR_PREC (default 512 bits).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd.h>
typedef dd_real T; typedef dd_complex CT; static const int NW = 2;
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td.h>
typedef td_real T; typedef td_complex CT; static const int NW = 3;
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include <mpblas_qd.h>
typedef qd_real T; typedef qd_complex CT; static const int NW = 4;
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
#include <mpblas_mpfr.h>
typedef mpfr_class T; typedef mpc_class CT;
#define IS_MPFR 1
#endif

static std::mt19937_64 g(12345);
static std::uniform_real_distribution<double> u(-1.0, 1.0);

#ifdef IS_MPFR
static T rnd() { // 8 doubles of 53 bits -> fills up to ~424+ bits; enough random low bits for timing
    T x = u(g); T s = 1.0;
    for (int i = 1; i < 10; i++) { s = s * T(0x1p-53); x = x + T(u(g)) * s; }
    return x;
}
static double todbl(const T &x) { return mpfr_get_d(x.get_mpfr_t(), MPFR_RNDN); }
static CT mkc(const T &a, const T &b) { return CT(a, b); }
static void rd(T &x, FILE *f) { std::vector<char> s(4096); fscanf(f, "%4095s", s.data()); mpfr_set_str(x.get_mpfr_t(), s.data(), 16, MPFR_RNDN); }
static void wr(const T &x, FILE *f) { mpfr_out_str(f, 16, 0, x.get_mpfr_t(), MPFR_RNDN); fputc('\n', f); }
static T re(const CT &z) { return T(z.real()); }
static T im(const CT &z) { return T(z.imag()); }
#else
static T rnd() {
    T x = u(g); double s = 1.0;
    for (int i = 1; i < NW; i++) { s *= 0x1p-53; x += T(u(g) * s); }
    return x;
}
static double todbl(const T &x) { double s = 0; for (int i = NW - 1; i >= 0; i--) s += x.x[i]; return s; }
static CT mkc(const T &a, const T &b) { return CT(a, b); }
static void rd(T &x, FILE *f) { fread(x.x, sizeof(double), NW, f); }
static void wr(const T &x, FILE *f) { fwrite(x.x, sizeof(double), NW, f); }
static T re(const CT &z) { return z.real(); }
static T im(const CT &z) { return z.imag(); }
#endif

int main(int argc, char **argv) {
#ifdef IS_MPFR
    mpfr_set_default_prec(getenv("MPFR_PREC") ? atol(getenv("MPFR_PREC")) : 512);
#endif
    if (argc >= 6 && !strcmp(argv[1], "cmp")) {
        bool cplx = argv[2][0] == 'C';
        long n = atol(argv[3]);
        size_t cnt = (size_t)n * n * (cplx ? 2 : 1);
        FILE *fa = fopen(argv[4], "rb"), *fb = fopen(argv[5], "rb");
        double mref = 0, md = 0;
        T a, b;
        for (size_t i = 0; i < cnt; i++) {
            rd(a, fa); rd(b, fb);
            double r = std::fabs(todbl(a)); if (r > mref) mref = r;
            T d = a - b; double dd = std::fabs(todbl(d)); if (dd > md) md = dd;
        }
        fclose(fa); fclose(fb);
        printf("max|diff|=%.3e max|ref|=%.3e rel=%.3e\n", md, mref, md / mref);
        return 0;
    }
    bool cplx = argv[1][0] == 'C';
    mplapackint n = atol(argv[2]);
    int reps = argc > 3 ? atoi(argv[3]) : 3;
    const char *dump = argc > 4 ? argv[4] : nullptr;
    double best = 1e300;
    if (!cplx) {
        std::vector<T> A(n * n), B(n * n), C(n * n);
        for (auto &x : A) x = rnd();
        for (auto &x : B) x = rnd();
        for (int r = 0; r < reps; r++) {
            auto t0 = std::chrono::steady_clock::now();
            Rgemm("N", "N", n, n, n, T(1.0), A.data(), n, B.data(), n, T(0.0), C.data(), n);
            double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (t < best) best = t;
        }
        if (dump) { FILE *f = fopen(dump, "wb"); for (auto &x : C) wr(x, f); fclose(f); }
        printf("Rgemm n=%ld best=%.4f s %.3f GFLOPS\n", (long)n, best, 2.0 * n * n * n / best * 1e-9);
    } else {
        std::vector<CT> A(n * n), B(n * n), C(n * n);
        for (auto &x : A) { T a = rnd(), b = rnd(); x = mkc(a, b); }
        for (auto &x : B) { T a = rnd(), b = rnd(); x = mkc(a, b); }
        CT one = mkc(T(1.0), T(0.0)), zero = mkc(T(0.0), T(0.0));
        for (int r = 0; r < reps; r++) {
            auto t0 = std::chrono::steady_clock::now();
            Cgemm("N", "N", n, n, n, one, A.data(), n, B.data(), n, zero, C.data(), n);
            double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (t < best) best = t;
        }
        if (dump) { FILE *f = fopen(dump, "wb"); for (auto &x : C) { wr(re(x), f); wr(im(x), f); } fclose(f); }
        printf("Cgemm n=%ld best=%.4f s %.3f GFLOPS\n", (long)n, best, 8.0 * n * n * n / best * 1e-9);
    }
    return 0;
}
