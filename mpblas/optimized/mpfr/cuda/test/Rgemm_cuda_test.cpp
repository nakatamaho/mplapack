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
 * Checks the fixed-precision CUDA Rgemm path against the CPU Rgemm of
 * libmplapack_mpfr_opt.  Linked with Rgemm_bridge_cuda.cpp and either
 * Rgemm_host_cuda.cpp (host emulation of the kernels, no GPU needed) or
 * Rgemm_device_cuda.cu (real GPU).
 *
 * For every eligible call the GPU path must return true and produce the
 * same MPFR value, limb for limb, as the CPU.  With the fixed-precision
 * kernels (512/1024 bits) a zero may differ in sign; the runtime-precision
 * kernels (other precisions, or all with MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1)
 * must match exactly, signed zeros included.
 * Ineligible calls must return false and leave C unchanged.
 */

#include <mpblas_mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace mplapack_mpfr_cuda {
int device_count();
}

bool Rgemm_mpfr_cuda(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc);

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
bool signed_zero_may_differ(int prec) { return !force_runtime && (prec == 512 || prec == 1024); }

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

void check_eligible(int prec, const char *ta, const char *tb, mplapackint m, mplapackint n, mplapackint k, double beta_d, int alpha_kind)
{
    bool nota = Mlsame_mpfr(ta, "N"), notb = Mlsame_mpfr(tb, "N");
    mplapackint nrowa = nota ? m : k, ncola = nota ? k : m;
    mplapackint nrowb = notb ? k : n, ncolb = notb ? n : k;
    mplapackint lda = nrowa + 3, ldb = nrowb + 1, ldc = m + 2;
    std::vector<mpfr_class> A((size_t)lda * ncola), B((size_t)ldb * ncolb), C0((size_t)ldc * n);
    fill(A);
    fill(B);
    fill(C0);
    mpfr_class alpha, beta = beta_d;
    if (alpha_kind == 0)
        alpha = 1.0;
    else if (alpha_kind == 1)
        alpha = -0.5;
    else
        do // Rgemm returns before the GPU path when alpha == 0
            random_value(alpha);
        while (alpha == 0.0);

    std::vector<mpfr_class> Cref(C0), Cgpu(C0);
    Rgemm(ta, tb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, Cref.data(), ldc);
    bool used = Rgemm_mpfr_cuda(nota, notb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, Cgpu.data(), ldc);
    if (!used) {
        std::printf("FAIL prec=%d %s%s m=%ld n=%ld k=%ld beta=%g alpha_kind=%d: GPU path not taken\n", prec, ta, tb, (long)m, (long)n, (long)k, beta_d, alpha_kind);
        failures++;
        return;
    }
    long bad = 0;
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < ldc; i++) {
            size_t t = i + j * ldc;
            checked++;
            if (!same_bits(Cref[t], Cgpu[t], signed_zero_may_differ(prec))) {
                if (bad == 0)
                    mpfr_printf("FAIL prec=%d %s%s m=%ld n=%ld k=%ld beta=%g alpha_kind=%d at (%ld,%ld): cpu=%.20Re gpu=%.20Re\n", prec, ta, tb, (long)m, (long)n, (long)k, beta_d, alpha_kind, (long)i, (long)j, Cref[t].get_mpfr_t(), Cgpu[t].get_mpfr_t());
                bad++;
            }
        }
    if (bad) {
        std::printf("  %ld mismatching elements\n", bad);
        failures++;
    }
}

// Calls that must not be taken by the GPU path; C must stay unchanged.
void check_ineligible(int prec)
{
    const mplapackint m = 20, n = 20, k = 20;
    std::vector<mpfr_class> A(m * k), B(k * n), C(m * n);
    fill(A);
    fill(B);
    fill(C);
    mpfr_class alpha = 1.5, beta = 0.25;
    struct {
        const char *what;
        int kind;
    } cases[] = {{"NaN in A", 0}, {"Inf in B", 1}, {"mixed precision in C", 2}, {"exponent near emax", 3}, {"rounding mode not RNDN", 4}};
    for (auto &c : cases) {
        std::vector<mpfr_class> A2(A), B2(B), C2(C);
        if (c.kind == 0)
            mpfr_set_nan(A2[5].get_mpfr_t());
        if (c.kind == 1)
            mpfr_set_inf(B2[7].get_mpfr_t(), -1);
        if (c.kind == 2)
            mpfr_set_prec(C2[3].get_mpfr_t(), prec + 64);
        if (c.kind == 3)
            mpfr_set_ui_2exp(A2[2].get_mpfr_t(), 1, mpfr_get_emax() - 1, MPFR_RNDN);
        std::vector<mpfr_class> before(C2);
        if (c.kind == 4)
            mpfr_set_default_rounding_mode(MPFR_RNDZ);
        bool used = Rgemm_mpfr_cuda(true, true, m, n, k, alpha, A2.data(), m, B2.data(), k, beta, C2.data(), m);
        mpfr_set_default_rounding_mode(MPFR_RNDN);
        bool unchanged = true;
        for (size_t t = 0; t < C2.size(); t++)
            if (mpfr_get_prec(C2[t].get_mpfr_t()) != mpfr_get_prec(before[t].get_mpfr_t()) || (!mpfr_nan_p(before[t].get_mpfr_t()) && !same_bits(C2[t], before[t])))
                unchanged = false;
        if (used || !unchanged) {
            std::printf("FAIL prec=%d ineligible case '%s': used=%d unchanged=%d\n", prec, c.what, (int)used, (int)unchanged);
            failures++;
        }
    }
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
    gmp_randseed_ui(rng, 20261002);

    const char *fr = std::getenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME");
    force_runtime = fr != NULL && fr[0] != '\0' && fr[0] != '0';

    // 512 and 1024 use the fixed-precision kernels (unless forced to the
    // runtime-precision ones); the others use the runtime-precision kernels.
    const int precs[] = {512, 1024, 64, 200, 256, 333, 2048};
    const char *trans[] = {"N", "T"};
    const double betas[] = {0.0, 1.0, 0.75, -1.25};
    const mplapackint shapes[][3] = {{1, 1, 1}, {7, 5, 3}, {37, 29, 41}, {16, 64, 8}};
    for (int prec : precs) {
        mpfrxx::set_default_precision_bits(prec);
        for (const char *ta : trans)
            for (const char *tb : trans)
                for (double beta : betas)
                    for (int alpha_kind = 0; alpha_kind < 3; alpha_kind++)
                        for (auto &sh : shapes)
                            check_eligible(prec, ta, tb, sh[0], sh[1], sh[2], beta, alpha_kind);
        check_ineligible(prec);
    }

    // Operands of a precision other than the default one: never eligible.
    mpfrxx::set_default_precision_bits(256);
    {
        std::vector<mpfr_class> A(64), B(64), C(64);
        fill(A);
        fill(B);
        fill(C);
        mpfr_class alpha = 1.0, beta = 1.0;
        mpfrxx::set_default_precision_bits(320);
        bool used = Rgemm_mpfr_cuda(true, true, 8, 8, 8, alpha, A.data(), 8, B.data(), 8, beta, C.data(), 8);
        if (used) {
            std::printf("FAIL operands at 256 bits with default precision 320 were taken by the GPU path\n");
            failures++;
        }
    }

    gmp_randclear(rng);
    std::printf("%s: %ld elements compared, %ld zeros differing only in sign, %d failures\n", failures ? "FAILED" : "PASSED", checked, signed_zero_diffs, failures);
    return failures ? 1 : 0;
}
