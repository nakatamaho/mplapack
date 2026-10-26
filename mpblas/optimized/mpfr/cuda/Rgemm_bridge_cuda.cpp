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
 * Bridge between mpfr_class matrices and the CUDA Rgemm.
 *
 * The entry point Rgemm_mpfr_cuda() (include/mpblas_mpfr_cuda.h) returns
 * true when it has computed C on the GPU, and false (leaving C untouched)
 * when the call is not eligible, so that Rgemm falls back to the CPU code.
 * A call is eligible when
 *   - the MPFR default precision, alpha, beta and every referenced element of
 *     A, B and C have the same precision p,
 *   - the MPFR default rounding mode is MPFR_RNDN,
 *   - every referenced value is finite (C is not read when beta == 0),
 *   - exponents are small enough that no intermediate can leave the MPFR
 *     exponent range, and every result fits in it,
 *   - m*n*k >= MPLAPACK_MPFR_CUDA_MIN_MNK (default 32768), and
 *   - MPLAPACK_MPFR_CUDA is not set to 0.
 *
 * p = 256, 512, 768, 1024 or 2048 uses the fixed-precision kernels (cu_fp::cu_freal<p>,
 * Rgemm_device_cuda.cu); the result equals the CPU Rgemm except that a zero
 * is always +0.  Any other p uses the runtime-precision kernels (cu_mpfr,
 * Rgemm_rt_cuda.cu), whose result equals the CPU Rgemm exactly.
 * MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1 sends those calls to the
 * runtime-precision kernels too; MPLAPACK_MPFR_CUDA_RUNTIME=0 disables them.
 *
 * MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF=N (N > 0) computes fixed-precision calls
 * with min(m, n, k) > N by the Winograd variant of Strassen's algorithm
 * (Rgemm_winograd_cuda.h), recursing down to blocks of size N.  The result
 * then differs from the CPU Rgemm within a normwise error bound.  It is off
 * by default.
 */

#include <mpblas_mpfr_cuda.h>
#include "mpfr_pack_cuda.h"
#include "../../common/Rgemm_entry_common.h"
#include "Rgemm_winograd_cuda.h"

namespace {

using namespace mplapack_mpfr_cuda_host;

template <int PB> bool run(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc)
{
    typedef real_t<PB> F;
    const mpfr_exp_t bound = input_exponent_bound();
    gemm_shape s;
    s.nota = nota;
    s.notb = notb;
    s.m = m;
    s.n = n;
    s.k = k;
    s.lda = nota ? m : k;
    s.ldb = notb ? k : n;
    s.ldc = m;
    s.beta_zero = (beta == 0.0);
    s.beta_one = (beta == 1.0);

    F a, b;
    if (!pack<PB>(alpha.get_mpfr_t(), a, bound) || !pack<PB>(beta.get_mpfr_t(), b, bound))
        return false;
    const long cutoff = winograd_cutoff();
    const bool use_winograd = cutoff > 0 && m > cutoff && n > cutoff && k > cutoff;
    std::vector<F> pA, pB, pC;
    if (use_winograd) {
        // op(A) and op(B) explicitly: m x k and k x n
        if (!pack_op_matrix<PB>(A, lda, !nota, m, k, pA, bound) || !pack_op_matrix<PB>(B, ldb, !notb, k, n, pB, bound))
            return false;
    } else if (!pack_matrix<PB>(A, lda, s.lda, nota ? k : m, pA, bound) || !pack_matrix<PB>(B, ldb, s.ldb, notb ? n : k, pB, bound)) {
        return false;
    }
    if (s.beta_zero) {
        // C is overwritten without being read; only its precision matters.
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                if (mpfr_get_prec(C[i + j * ldc].get_mpfr_t()) != PB)
                    return false;
        pC.resize((size_t)m * n);
    } else if (!pack_matrix<PB>(C, ldc, m, n, pC, bound)) {
        return false;
    }

    const int rc = use_winograd ? mplapack_mpfr_cuda::gemm_winograd<PB>(s, a, b, pA.data(), pB.data(), pC.data(), cutoff) : mplapack_mpfr_cuda::gemm<PB>(s, a, b, pA.data(), pB.data(), pC.data());
    if (rc != 0)
        return false;
    for (size_t t = 0; t < pC.size(); t++)
        if (!fits<PB>(pC[t]))
            return false;
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            unpack<PB>(pC[i + j * m], C[i + j * ldc].get_mpfr_t());
    return true;
}

bool run_rt(mpfr_prec_t prec, bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc)
{
    const int nl = (int)((prec + 63) / 64);
    const mpfr_exp_t bound = rt_input_exponent_bound();
    gemm_shape s;
    s.nota = nota;
    s.notb = notb;
    s.m = m;
    s.n = n;
    s.k = k;
    s.lda = nota ? m : k;
    s.ldb = notb ? k : n;
    s.ldc = m;
    s.beta_zero = (beta == 0.0);
    s.beta_one = (beta == 1.0);

    rt_storage sa, sb, pA, pB, pC;
    sa.resize(1, nl);
    sb.resize(1, nl);
    if (!rt_pack_value(alpha.get_mpfr_t(), prec, nl, sa, 0, bound) || !rt_pack_value(beta.get_mpfr_t(), prec, nl, sb, 0, bound))
        return false;
    const long ncola = nota ? k : m, ncolb = notb ? n : k;
    if (!rt_pack_matrix(A, lda, s.lda, ncola, prec, nl, pA, bound) || !rt_pack_matrix(B, ldb, s.ldb, ncolb, prec, nl, pB, bound))
        return false;
    if (s.beta_zero) {
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                if (mpfr_get_prec(C[i + j * ldc].get_mpfr_t()) != prec)
                    return false;
        pC.resize((size_t)m * n, nl);
    } else if (!rt_pack_matrix(C, ldc, m, n, prec, nl, pC, bound)) {
        return false;
    }

    if (mplapack_mpfr_cuda::gemm_rt(s, (long)prec, sa.pack, sb.pack, pA.pack, (size_t)s.lda * ncola, pB.pack, (size_t)s.ldb * ncolb, pC.pack) != 0)
        return false;
    for (size_t t = 0; t < (size_t)m * n; t++)
        if (!rt_fits(pC, t))
            return false;
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            rt_unpack(pC, (size_t)(i + j * m), prec, nl, C[i + j * ldc].get_mpfr_t());
    return true;
}

} // namespace

// Entry point (include/mpblas_mpfr_cuda.h).  Invalid arguments, the quick
// returns of Rgemm and alpha == 0 are left to the caller.
bool Rgemm_mpfr_cuda(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *A, mplapackint const lda, mpfr_class *B, mplapackint const ldb, mpfr_class const beta, mpfr_class *C, mplapackint const ldc)
{
    bool nota, notb;
    if (!mplapack_gemm::gemm_entry_args(transa, transb, m, n, k, lda, ldb, ldc, alpha == 0.0, beta == 1.0, nota, notb) || alpha == 0.0)
        return false;
    if (!cuda_enabled() || (double)m * n * k < (double)min_mnk())
        return false;
    if (mpfr_get_default_rounding_mode() != MPFR_RNDN)
        return false;
    // The CPU code rounds temporaries to the default precision.
    const mpfr_prec_t prec = mpfrxx::default_precision_bits();
    if (!force_runtime()) {
        switch (prec) {
#define MPLAPACK_MPFR_CUDA_DISPATCH(PB) \
    case PB:                            \
        return run<PB>(nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_DISPATCH)
#undef MPLAPACK_MPFR_CUDA_DISPATCH
        default:
            break;
        }
    }
    if (!runtime_enabled())
        return false;
    return run_rt(prec, nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
