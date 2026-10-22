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
 * Conversion between mpfr_class values and the packed operands of the MPFR
 * CUDA BLAS (fixed precision: cu_fp::cu_freal<PB>; runtime precision:
 * rt_pack of limbs / exponents / kinds), and the settings shared by the
 * bridges (environment variables, exponent bounds).  Host code only.
 */

#ifndef _MPLAPACK_MPFR_PACK_CUDA_H_
#define _MPLAPACK_MPFR_PACK_CUDA_H_

#include <mpblas_mpfr.h>
#include <cstdlib>
#include <vector>
#include "Rgemm_kernel_cuda.h"
#include "Rgemm_rt_cuda.h"

static_assert(GMP_NUMB_BITS == 64 && sizeof(mp_limb_t) == sizeof(cu_fp::cu_limb), "mpfr cuda BLAS requires 64-bit GMP limbs");
static_assert(MPFR_ZERO_KIND == mplapack_mpfr_cuda::RT_ZERO_KIND && MPFR_REGULAR_KIND == mplapack_mpfr_cuda::RT_REGULAR_KIND, "MPFR custom kinds differ");

namespace mplapack_mpfr_cuda_host {

using mplapack_mpfr_cuda::gemm_shape;
using mplapack_mpfr_cuda::real_t;

inline long env_long(const char *name, long fallback)
{
    const char *v = std::getenv(name);
    if (v == NULL || v[0] == '\0')
        return fallback;
    char *end = NULL;
    long r = std::strtol(v, &end, 10);
    return (end != v && *end == '\0') ? r : fallback;
}

inline bool cuda_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA", 1) != 0;
    return enabled;
}

inline long min_mnk()
{
    static const long v = env_long("MPLAPACK_MPFR_CUDA_MIN_MNK", 32768);
    return v;
}

inline bool runtime_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA_RUNTIME", 1) != 0;
    return enabled;
}

inline long winograd_cutoff()
{
    static const long v = env_long("MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF", 0);
    return v;
}

inline bool force_runtime()
{
    static const bool forced = env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME", 0) != 0;
    return forced;
}

// Bound on |exponent| of the inputs: alpha*op(B)*op(A) multiplies three
// values and sums at most k terms, so 4*bound + 64 must stay inside the range.
inline mpfr_exp_t input_exponent_bound()
{
    mpfr_exp_t emax = mpfr_get_emax();
    mpfr_exp_t emin = -mpfr_get_emin();
    mpfr_exp_t e = emax < emin ? emax : emin;
    return (e - 128) / 4;
}

template <int PB> bool pack(mpfr_srcptr x, real_t<PB> &r, mpfr_exp_t bound)
{
    if (mpfr_get_prec(x) != PB)
        return false;
    if (mpfr_zero_p(x)) {
        r.set_zero();
        return true;
    }
    if (!mpfr_regular_p(x))
        return false;
    mpfr_exp_t e = mpfr_get_exp(x);
    if (e > bound || e < -bound)
        return false;
    r.sign = mpfr_signbit(x) ? -1 : 1;
    r.exp = (long)e;
    const mp_limb_t *d = static_cast<const mp_limb_t *>(mpfr_custom_get_significand(x));
    for (int i = 0; i < real_t<PB>::N; i++)
        r.m[i] = d[i];
    return true;
}

template <int PB> bool fits(const real_t<PB> &r) { return r.is_zero() || (r.exp >= mpfr_get_emin() && r.exp <= mpfr_get_emax()); }

template <int PB> void unpack(const real_t<PB> &r, mpfr_ptr x)
{
    if (r.is_zero()) {
        mpfr_set_zero(x, 1);
        return;
    }
    mp_limb_t limbs[real_t<PB>::N];
    for (int i = 0; i < real_t<PB>::N; i++)
        limbs[i] = r.m[i];
    mpfr_t v;
    mpfr_custom_init_set(v, r.sign > 0 ? MPFR_REGULAR_KIND : -MPFR_REGULAR_KIND, (mpfr_exp_t)r.exp, PB, limbs);
    mpfr_set(x, v, MPFR_RNDN); // exact: x has precision PB
}

// Packs rows x cols elements of X (leading dimension ld) into out (leading dimension rows).
template <int PB> bool pack_matrix(const mpfr_class *X, mplapackint ld, long rows, long cols, std::vector<real_t<PB>> &out, mpfr_exp_t bound)
{
    out.resize((size_t)rows * cols);
    for (long j = 0; j < cols; j++)
        for (long i = 0; i < rows; i++)
            if (!pack<PB>(X[i + j * ld].get_mpfr_t(), out[i + j * rows], bound))
                return false;
    return true;
}

// Packs op(X) (rows x cols) with leading dimension rows; trans: X is stored transposed.
template <int PB> bool pack_op_matrix(const mpfr_class *X, mplapackint ld, bool trans, long rows, long cols, std::vector<real_t<PB>> &out, mpfr_exp_t bound)
{
    out.resize((size_t)rows * cols);
    for (long j = 0; j < cols; j++)
        for (long i = 0; i < rows; i++)
            if (!pack<PB>((trans ? X[j + i * ld] : X[i + j * ld]).get_mpfr_t(), out[i + j * rows], bound))
                return false;
    return true;
}

// ---- runtime precision (cu_mpfr) ----

// The device MPFR uses its default exponent range, which may be narrower
// than the host one.
inline mpfr_exp_t rt_input_exponent_bound()
{
    mpfr_exp_t e = input_exponent_bound();
    mpfr_exp_t d = (mplapack_mpfr_cuda::RT_DEVICE_EMAX - 128) / 4;
    return e < d ? e : d;
}

struct rt_storage {
    std::vector<unsigned long long> limbs;
    std::vector<long> exp;
    std::vector<signed char> kind;
    mplapack_mpfr_cuda::rt_pack pack;
    void resize(size_t count, int nl)
    {
        limbs.assign(count * nl + 1, 0);
        exp.assign(count + 1, 0);
        kind.assign(count + 1, 0);
        pack.limbs = limbs.data();
        pack.exp = exp.data();
        pack.kind = kind.data();
    }
};

inline bool rt_pack_value(mpfr_srcptr x, mpfr_prec_t prec, int nl, rt_storage &st, size_t idx, mpfr_exp_t bound)
{
    if (mpfr_get_prec(x) != prec)
        return false;
    int kind = mpfr_custom_get_kind(x);
    if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND) {
        st.kind[idx] = (signed char)kind;
        return true;
    }
    if (kind != MPFR_REGULAR_KIND && kind != -MPFR_REGULAR_KIND)
        return false;
    mpfr_exp_t e = mpfr_get_exp(x);
    if (e > bound || e < -bound)
        return false;
    st.kind[idx] = (signed char)kind;
    st.exp[idx] = (long)e;
    const mp_limb_t *d = static_cast<const mp_limb_t *>(mpfr_custom_get_significand(x));
    for (int t = 0; t < nl; t++)
        st.limbs[idx * nl + t] = d[t];
    return true;
}

inline bool rt_pack_matrix(const mpfr_class *X, mplapackint ld, long rows, long cols, mpfr_prec_t prec, int nl, rt_storage &st, mpfr_exp_t bound)
{
    st.resize((size_t)rows * cols, nl);
    for (long j = 0; j < cols; j++)
        for (long i = 0; i < rows; i++)
            if (!rt_pack_value(X[i + j * ld].get_mpfr_t(), prec, nl, st, (size_t)(i + j * rows), bound))
                return false;
    return true;
}

inline bool rt_fits(const rt_storage &st, size_t idx)
{
    int kind = st.kind[idx];
    if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND)
        return true;
    if (kind != MPFR_REGULAR_KIND && kind != -MPFR_REGULAR_KIND)
        return false;
    return st.exp[idx] >= mpfr_get_emin() && st.exp[idx] <= mpfr_get_emax();
}

inline void rt_unpack(const rt_storage &st, size_t idx, mpfr_prec_t prec, int nl, mpfr_ptr x)
{
    int kind = st.kind[idx];
    if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND) {
        mpfr_set_zero(x, kind > 0 ? 1 : -1);
        return;
    }
    mpfr_t v;
    mpfr_custom_init_set(v, kind, (mpfr_exp_t)st.exp[idx], prec, const_cast<unsigned long long *>(&st.limbs[idx * nl]));
    mpfr_set(x, v, MPFR_RNDN); // exact: x has precision prec
}

} // namespace mplapack_mpfr_cuda_host

#endif
