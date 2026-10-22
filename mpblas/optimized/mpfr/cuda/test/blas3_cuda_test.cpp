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
 * Checks the MPFR Level 3 BLAS on CUDA (Rsyrk, Rsyr2k, Rsymm, Rtrmm, Rtrsm)
 * against the reference implementation of libmplapack_mpfr_opt, bit for bit,
 * for every side / uplo / trans / diag, alpha in {1, -1/2, random},
 * beta in {0, 1, 3/4, -5/4}, odd sizes and leading dimensions, at the
 * fixed precisions (cu_freal engine; zeros may differ in sign) and at other
 * precisions (cu_mpfr engine; exact).  Linked with blas3_bridge_cuda.cpp,
 * blas3_device_cuda.cu and Rgemm_rt_cuda.cu; with
 * MPLAPACK_MPFR_CUDA_HOST_EMULATION the kernels run on the host.
 * Elements outside the referenced part of the result must be untouched.
 */

#include <mpblas_mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace mplapack_mpfr_cuda {
int device_count();
}

bool Rsyrk_mpfr_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class &beta, mpfr_class *c, mplapackint ldc);
bool Rsyr2k_mpfr_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class *b, mplapackint ldb, const mpfr_class &beta, mpfr_class *c, mplapackint ldc);
bool Rsymm_mpfr_cuda(const char *side, const char *uplo, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class *b, mplapackint ldb, const mpfr_class &beta, mpfr_class *c, mplapackint ldc);
bool Rtrmm_mpfr_cuda(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, mpfr_class *b, mplapackint ldb);
bool Rtrsm_mpfr_cuda(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, mpfr_class *b, mplapackint ldb);

namespace {

gmp_randstate_t rng;
int failures = 0;
long checked = 0;
long signed_zero_diffs = 0;

// Random value of mixed magnitude: a uniform significand scaled by 2^e with
// e in [-40, 40]; about 1 in 8 values is an exact zero and about 1 in 8 is a
// small integer, so cancellations and exact sums occur.
void random_value(mpfr_class &x)
{
    unsigned long r = gmp_urandomm_ui(rng, 8);
    if (r == 0) {
        x = 0.0;
        return;
    }
    if (r == 1) {
        x = (long)gmp_urandomm_ui(rng, 7) - 3;
        return;
    }
    mpfr_urandomb(x.get_mpfr_t(), rng);
    long e = (long)gmp_urandomm_ui(rng, 81) - 40;
    mpfr_mul_2si(x.get_mpfr_t(), x.get_mpfr_t(), e, MPFR_RNDN);
    if (gmp_urandomm_ui(rng, 2))
        mpfr_neg(x.get_mpfr_t(), x.get_mpfr_t(), MPFR_RNDN);
}

bool force_runtime = false;

// The fixed-precision kernels have no signed zero.
bool signed_zero_may_differ(int prec) { return !force_runtime && (prec == 256 || prec == 512 || prec == 768 || prec == 1024 || prec == 2048); }

bool same_bits(const mpfr_class &a, const mpfr_class &b, bool zero_sign_may_differ = false)
{
    mpfr_srcptr x = a.get_mpfr_t(), y = b.get_mpfr_t();
    if (mpfr_get_prec(x) != mpfr_get_prec(y))
        return false;
    if (mpfr_zero_p(x) && mpfr_zero_p(y)) {
        if (mpfr_signbit(x) != mpfr_signbit(y)) {
            if (!zero_sign_may_differ)
                return false;
            signed_zero_diffs++;
        }
        return true;
    }
    if (mpfr_nan_p(x) || mpfr_nan_p(y))
        return mpfr_nan_p(x) && mpfr_nan_p(y);
    if (mpfr_inf_p(x) || mpfr_inf_p(y))
        return mpfr_inf_p(x) && mpfr_inf_p(y) && mpfr_signbit(x) == mpfr_signbit(y);
    if (!mpfr_regular_p(x) || !mpfr_regular_p(y))
        return false;
    if (mpfr_signbit(x) != mpfr_signbit(y) || mpfr_get_exp(x) != mpfr_get_exp(y))
        return false;
    size_t nlimbs = (mpfr_get_prec(x) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    return std::memcmp(mpfr_custom_get_significand(x), mpfr_custom_get_significand(y), nlimbs * sizeof(mp_limb_t)) == 0;
}


void fill(std::vector<mpfr_class> &v)
{
    for (size_t i = 0; i < v.size(); i++)
        random_value(v[i]);
}

// A triangular matrix that keeps Rtrsm well conditioned: diagonal in
// [1, 2) (sign random), off-diagonal entries scaled down by 1/n.
void fill_triangular(std::vector<mpfr_class> &A, long n, long lda)
{
    fill(A);
    for (long j = 0; j < n; j++)
        for (long i = 0; i < n; i++) {
            mpfr_class &x = A[i + j * lda];
            if (i == j) {
                mpfr_urandomb(x.get_mpfr_t(), rng);
                x += 1.0;
                if (gmp_urandomm_ui(rng, 2))
                    x = -x;
            } else {
                x /= (double)(2 * n);
            }
        }
}

mpfr_class pick_alpha(int kind)
{
    mpfr_class a;
    if (kind == 0)
        a = 1.0;
    else if (kind == 1)
        a = -0.5;
    else
        do
            random_value(a);
        while (a == 0.0);
    return a;
}

void compare(const std::string &what, int prec, const std::vector<mpfr_class> &ref, const std::vector<mpfr_class> &got, bool used)
{
    if (!used) {
        std::printf("FAIL prec=%d %s: GPU path not taken\n", prec, what.c_str());
        failures++;
        return;
    }
    long bad = 0;
    for (size_t t = 0; t < ref.size(); t++) {
        checked++;
        if (!same_bits(ref[t], got[t], signed_zero_may_differ(prec))) {
            if (bad == 0)
                mpfr_printf("FAIL prec=%d %s at %ld: cpu=%.20Re gpu=%.20Re\n", prec, what.c_str(), (long)t, ref[t].get_mpfr_t(), got[t].get_mpfr_t());
            bad++;
        }
    }
    if (bad) {
        std::printf("  %ld mismatching elements\n", bad);
        failures++;
    }
}

const char *yn(bool b, const char *y, const char *n) { return b ? y : n; }

void check_syrk(int prec, bool two, const char *uplo, const char *trans, long n, long k, double beta_d, int ak)
{
    const bool notr = trans[0] == 'N';
    const long nrowa = notr ? n : k, ncola = notr ? k : n, lda = nrowa + 2, ldc = n + 1;
    std::vector<mpfr_class> A(lda * ncola), B(lda * ncola), C(ldc * n);
    fill(A);
    fill(B);
    fill(C);
    mpfr_class alpha = pick_alpha(ak), beta = beta_d;
    std::vector<mpfr_class> R(C), G(C);
    bool used;
    if (two) {
        Rsyr2k(uplo, trans, n, k, alpha, A.data(), lda, B.data(), lda, beta, R.data(), ldc);
        used = Rsyr2k_mpfr_cuda(uplo, trans, n, k, alpha, A.data(), lda, B.data(), lda, beta, G.data(), ldc);
    } else {
        Rsyrk(uplo, trans, n, k, alpha, A.data(), lda, beta, R.data(), ldc);
        used = Rsyrk_mpfr_cuda(uplo, trans, n, k, alpha, A.data(), lda, beta, G.data(), ldc);
    }
    char buf[160];
    std::snprintf(buf, sizeof(buf), "%s %s%s n=%ld k=%ld beta=%g alpha_kind=%d", two ? "Rsyr2k" : "Rsyrk", uplo, trans, n, k, beta_d, ak);
    compare(buf, prec, R, G, used);
}

void check_symm(int prec, const char *side, const char *uplo, long m, long n, double beta_d, int ak)
{
    const long na = side[0] == 'L' ? m : n, lda = na + 1, ldb = m + 2, ldc = m + 3;
    std::vector<mpfr_class> A(lda * na), B(ldb * n), C(ldc * n);
    fill(A);
    fill(B);
    fill(C);
    mpfr_class alpha = pick_alpha(ak), beta = beta_d;
    std::vector<mpfr_class> R(C), G(C);
    Rsymm(side, uplo, m, n, alpha, A.data(), lda, B.data(), ldb, beta, R.data(), ldc);
    bool used = Rsymm_mpfr_cuda(side, uplo, m, n, alpha, A.data(), lda, B.data(), ldb, beta, G.data(), ldc);
    char buf[160];
    std::snprintf(buf, sizeof(buf), "Rsymm %s%s m=%ld n=%ld beta=%g alpha_kind=%d", side, uplo, m, n, beta_d, ak);
    compare(buf, prec, R, G, used);
}

void check_trxm(int prec, bool solve, const char *side, const char *uplo, const char *transa, const char *diag, long m, long n, int ak)
{
    const long na = side[0] == 'L' ? m : n, lda = na + 2, ldb = m + 1;
    std::vector<mpfr_class> A(lda * na), B(ldb * n);
    if (solve)
        fill_triangular(A, na, lda);
    else
        fill(A);
    fill(B);
    mpfr_class alpha = pick_alpha(ak);
    std::vector<mpfr_class> R(B), G(B);
    bool used;
    if (solve) {
        Rtrsm(side, uplo, transa, diag, m, n, alpha, A.data(), lda, R.data(), ldb);
        used = Rtrsm_mpfr_cuda(side, uplo, transa, diag, m, n, alpha, A.data(), lda, G.data(), ldb);
    } else {
        Rtrmm(side, uplo, transa, diag, m, n, alpha, A.data(), lda, R.data(), ldb);
        used = Rtrmm_mpfr_cuda(side, uplo, transa, diag, m, n, alpha, A.data(), lda, G.data(), ldb);
    }
    char buf[160];
    std::snprintf(buf, sizeof(buf), "%s %s%s%s%s m=%ld n=%ld alpha_kind=%d", solve ? "Rtrsm" : "Rtrmm", side, uplo, transa, diag, m, n, ak);
    // Rtrsm runs on the cu_mpfr engine at every precision: exact.
    const bool saved = force_runtime;
    if (solve)
        force_runtime = true;
    compare(buf, prec, R, G, used);
    force_runtime = saved;
}

} // namespace

int main()
{
    if (mplapack_mpfr_cuda::device_count() == 0) {
        std::printf("SKIPPED: no CUDA device\n");
        return 77;
    }
    setenv("MPLAPACK_MPFR_CUDA_MIN_MNK", "0", 1);
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 20261022);
    const char *fr = std::getenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME");
    force_runtime = fr != NULL && fr[0] != '\0' && fr[0] != '0';

    const int precs[] = {256, 512, 1024, 2048, 768, 64, 200, 333};
    const double betas[] = {0.0, 1.0, 0.75, -1.25};
    const char *ul[] = {"U", "L"}, *tr[] = {"N", "T"}, *sd[] = {"L", "R"}, *dg[] = {"N", "U"};
    const long sizes[][2] = {{1, 1}, {7, 5}, {13, 17}};
    for (int prec : precs) {
        mpfrxx::set_default_precision_bits(prec);
        for (auto &z : sizes)
            for (int ak = 0; ak < 3; ak++) {
                for (const char *u : ul)
                    for (const char *t : tr)
                        for (double b : betas) {
                            check_syrk(prec, false, u, t, z[0], z[1], b, ak);
                            check_syrk(prec, true, u, t, z[0], z[1], b, ak);
                        }
                for (const char *s : sd)
                    for (const char *u : ul) {
                        for (double b : betas)
                            check_symm(prec, s, u, z[0], z[1], b, ak);
                        for (const char *t : tr)
                            for (const char *d : dg) {
                                check_trxm(prec, false, s, u, t, d, z[0], z[1], ak);
                                check_trxm(prec, true, s, u, t, d, z[0], z[1], ak);
                            }
                    }
            }
    }
    gmp_randclear(rng);
    std::printf("%s: %ld elements compared, %ld zeros differing only in sign, %d failures\n", failures ? "FAILED" : "PASSED", checked, signed_zero_diffs, failures);
    return failures ? 1 : 0;
}
