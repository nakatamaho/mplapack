// MPFR Rgemm("N","N",n,n,n) sweep at a fixed precision: n = powers of two and
// multiples of 13 up to NMAX.  A, B (NMAX x NMAX, uniform in [-1,1) at full
// precision) are generated once; each n uses their leading n x n blocks
// (lda = ldb = NMAX).  Best of REPS for n <= 2000, of REPS_LARGE above;
// sizes below NSTART are skipped (to resume).  Output: n,seconds,gflops
// usage: mpfr_sweep [PREC=512] [NMAX=2000] [REPS=3] [REPS_LARGE=REPS] [NSTART=1]
#include <mpblas_mpfr.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char **argv) {
    const long prec = argc > 1 ? atol(argv[1]) : 512;
    const long NMAX = argc > 2 ? atol(argv[2]) : 2000;
    const int reps = argc > 3 ? atoi(argv[3]) : 3;
    const int reps_large = argc > 4 ? atoi(argv[4]) : reps;
    const long nstart = argc > 5 ? atol(argv[5]) : 1;
    mpfr_set_default_prec(prec);
    std::vector<long> sizes;
    for (long n = 1; n <= NMAX; n *= 2) sizes.push_back(n);
    for (long n = 13; n <= NMAX; n += 13) sizes.push_back(n);
    std::sort(sizes.begin(), sizes.end());
    sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());

    auto t0 = std::chrono::steady_clock::now();
    std::vector<mpfr_class> A((size_t)NMAX * NMAX), B((size_t)NMAX * NMAX), C((size_t)NMAX * NMAX);
    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 12345);
    for (auto *M : {&A, &B})
        for (auto &x : *M) {
            mpfr_urandomb(x.get_mpfr_t(), rs);          // [0,1)
            mpfr_mul_2ui(x.get_mpfr_t(), x.get_mpfr_t(), 1, MPFR_RNDN);
            mpfr_sub_ui(x.get_mpfr_t(), x.get_mpfr_t(), 1, MPFR_RNDN); // [-1,1)
        }
    fprintf(stderr, "generated %ldx%ld in %.1f s\n", NMAX, NMAX, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    const mpfr_class one = 1.0, zero = 0.0;
    if (nstart <= 1)
        printf("n,seconds,gflops\n");
    for (long n : sizes) {
        if (n < nstart)
            continue;
        double best = 1e300;
        const int nrep = n <= 2000 ? reps : reps_large;
        for (int r = 0; r < nrep; r++) {
            auto t1 = std::chrono::steady_clock::now();
            Rgemm("N", "N", n, n, n, one, A.data(), NMAX, B.data(), NMAX, zero, C.data(), NMAX);
            double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
            best = std::min(best, t);
        }
        printf("%ld,%.6g,%.6g\n", n, best, 2.0 * n * n * n / best * 1e-9);
        fflush(stdout);
    }
    return 0;
}
