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
 * Checks the Winograd (Strassen) path of the CUDA Rgemm
 * (MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF), at 256, 512, 768, 1024 and 2048 bits:
 *
 * 1. Integer matrices: no rounding occurs, so Rgemm_mpfr_cuda must equal
 *    the CPU Rgemm exactly, for odd and even sizes, all transpose
 *    combinations, general alpha/beta and several recursion cutoffs.  This
 *    checks padding, block addressing, op() packing and the final
 *    alpha*P + beta*C.
 * 2. Random real matrices: the result must agree with the CPU Rgemm within
 *    a normwise bound, 2^(20-p) * k * max|alpha| * max|A| * max|B| + ...
 * 3. gemm_winograd (GPU, or host build) must equal the host reference
 *    backend bit for bit.
 */

#include <mpblas_mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include "Rgemm_winograd_cuda.h"

namespace mplapack_mpfr_cuda {
int device_count();
}

bool Rgemm_mpfr_cuda(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc);

namespace {

using mplapack_mpfr_cuda::gemm_shape;
using mplapack_mpfr_cuda::real_t;

gmp_randstate_t rng;
int failures = 0;
long checked = 0;

void random_int(mpfr_class &x) { x = (long)gmp_urandomm_ui(rng, 7) - 3; }

void random_real(mpfr_class &x)
{
    mpfr_urandomb(x.get_mpfr_t(), rng);
    mpfr_mul_2si(x.get_mpfr_t(), x.get_mpfr_t(), (long)gmp_urandomm_ui(rng, 9) - 4, MPFR_RNDN);
    if (gmp_urandomm_ui(rng, 2))
        mpfr_neg(x.get_mpfr_t(), x.get_mpfr_t(), MPFR_RNDN);
}

mpfr_class max_abs(const std::vector<mpfr_class> &v)
{
    mpfr_class r = 0.0;
    for (size_t i = 0; i < v.size(); i++)
        if (abs(v[i]) > r)
            r = abs(v[i]);
    return r;
}

// Rgemm on the CPU and through the GPU path; exact or bounded comparison.
void check(int prec, const char *ta, const char *tb, mplapackint m, mplapackint n, mplapackint k, bool integers)
{
    bool nota = Mlsame_mpfr(ta, "N"), notb = Mlsame_mpfr(tb, "N");
    mplapackint nrowa = nota ? m : k, ncola = nota ? k : m, nrowb = notb ? k : n, ncolb = notb ? n : k;
    mplapackint lda = nrowa + 2, ldb = nrowb + 1, ldc = m + 3;
    std::vector<mpfr_class> A((size_t)lda * ncola), B((size_t)ldb * ncolb), C0((size_t)ldc * n);
    for (auto &x : A)
        integers ? random_int(x) : random_real(x);
    for (auto &x : B)
        integers ? random_int(x) : random_real(x);
    for (auto &x : C0)
        integers ? random_int(x) : random_real(x);
    mpfr_class alpha = integers ? mpfr_class(2) : mpfr_class(0.75), beta = integers ? mpfr_class(-3) : mpfr_class(-1.25);

    std::vector<mpfr_class> Cref(C0), Cgpu(C0);
    Rgemm(ta, tb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, Cref.data(), ldc);
    if (!Rgemm_mpfr_cuda(nota, notb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, Cgpu.data(), ldc)) {
        std::printf("FAIL prec=%d %s%s %ldx%ldx%ld: GPU path not taken\n", prec, ta, tb, (long)m, (long)n, (long)k);
        failures++;
        return;
    }
    // bound for real data: 2^(20-p) * (k*|alpha|*max|A|*max|B| + |beta|*max|C|)
    mpfr_class bound = 0.0;
    if (!integers) {
        bound = (mpfr_class((double)k) * abs(alpha) * max_abs(A) * max_abs(B) + abs(beta) * max_abs(C0));
        mpfr_mul_2si(bound.get_mpfr_t(), bound.get_mpfr_t(), 20 - prec, MPFR_RNDN);
    }
    long bad = 0;
    for (size_t t = 0; t < Cref.size(); t++) {
        checked++;
        mpfr_class d = abs(Cref[t] - Cgpu[t]);
        if (integers ? (Cref[t] != Cgpu[t]) : (d > bound)) {
            if (bad == 0)
                mpfr_printf("FAIL prec=%d %s%s %ldx%ldx%ld %s: element %ld cpu=%.20Re gpu=%.20Re\n", prec, ta, tb, (long)m, (long)n, (long)k, integers ? "integers" : "reals", (long)t, Cref[t].get_mpfr_t(), Cgpu[t].get_mpfr_t());
            bad++;
        }
    }
    if (bad) {
        std::printf("  %ld mismatching elements\n", bad);
        failures++;
    }
}

// gemm_winograd (device or host build) against the host reference backend.
template <int PB> void check_backend(long m, long n, long k, long cutoff)
{
    typedef real_t<PB> F;
    std::vector<F> A((size_t)m * k), B((size_t)k * n), C((size_t)m * n), Cref;
    for (auto &x : A)
        x = F::from_double((double)gmp_urandomm_ui(rng, 1000001) / 1000.0 - 500.0);
    for (auto &x : B)
        x = F::from_double((double)gmp_urandomm_ui(rng, 1000001) / 3000.0 - 160.0);
    for (auto &x : C)
        x = F::from_double((double)gmp_urandomm_ui(rng, 1001) - 500.0);
    gemm_shape s;
    s.nota = s.notb = 1;
    s.m = m, s.n = n, s.k = k, s.lda = m, s.ldb = k, s.ldc = m;
    s.beta_zero = 0, s.beta_one = 0;
    F alpha = F::from_double(1.0 / 3.0), beta = F::from_double(-0.5);
    Cref = C;
    mplapack_mpfr_cuda::gemm_winograd_host<PB>(s, alpha, beta, A.data(), B.data(), Cref.data(), cutoff);
    if (mplapack_mpfr_cuda::gemm_winograd<PB>(s, alpha, beta, A.data(), B.data(), C.data(), cutoff) != 0) {
        std::printf("FAIL backend PB=%d %ldx%ldx%ld cutoff=%ld: gemm_winograd failed\n", PB, m, n, k, cutoff);
        failures++;
        return;
    }
    for (size_t t = 0; t < C.size(); t++) {
        checked++;
        if (cu_fp::cu_cmp<PB>(C[t], Cref[t]) != 0 || C[t].sign != Cref[t].sign) {
            std::printf("FAIL backend PB=%d %ldx%ldx%ld cutoff=%ld: element %ld differs\n", PB, m, n, k, cutoff, (long)t);
            failures++;
            return;
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
    setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME", "0", 1);
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 20261003);

    const char *trans[] = {"N", "T"};
    const mplapackint shapes[][3] = {{7, 5, 9}, {16, 16, 16}, {33, 17, 40}, {37, 29, 41}, {64, 64, 64}};
    const char *cutoffs[] = {"1", "2", "3", "8"};
    // The cutoff is read once per process, so each value runs in a child.
    const char *cut = std::getenv("MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF");
    if (cut == NULL) {
        int rc = 0;
        for (const char *c : cutoffs) {
            setenv("MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF", c, 1);
            char self[4000], cmd[4096];
            ssize_t len = readlink("/proc/self/exe", self, sizeof(self) - 1);
            if (len <= 0)
                return 1;
            self[len] = '\0';
            std::snprintf(cmd, sizeof(cmd), "'%s'", self);
            std::printf("--- cutoff %s\n", c);
            std::fflush(stdout);
            int r = std::system(cmd);
            if (r != 0)
                rc = 1;
        }
        return rc;
    }

    for (int prec : {256, 512, 768, 1024, 2048}) {
        mpfrxx::set_default_precision_bits(prec);
        for (const char *ta : trans)
            for (const char *tb : trans)
                for (auto &sh : shapes) {
                    check(prec, ta, tb, sh[0], sh[1], sh[2], true);
                    check(prec, ta, tb, sh[0], sh[1], sh[2], false);
                }
    }
    long cutoff = std::atol(cut);
    check_backend<256>(37, 29, 41, cutoff);
    check_backend<512>(37, 29, 41, cutoff);
    check_backend<768>(33, 17, 40, cutoff);
    check_backend<1024>(64, 64, 64, cutoff);
    check_backend<2048>(16, 16, 16, cutoff);

    gmp_randclear(rng);
    std::printf("%s: cutoff %ld, %ld elements compared, %d failures\n", failures ? "FAILED" : "PASSED", cutoff, checked, failures);
    return failures ? 1 : 0;
}
