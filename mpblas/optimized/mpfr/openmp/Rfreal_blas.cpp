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

// The entries of the fixed-precision MPFR BLAS (Rfreal_blas.h): the checks
// of the arithmetic context and of the CPU, then the driver of the precision
// (Rfreal_blas_<PB>.cpp).  Also the state of Rfreal_state.h.

#include <cstdlib>
#include <mpblas_mpfr.h>
#include "Rgemm_freal.h"
#include "Rfreal_blas.h"
#include "Rfreal_state.h"

thread_local mplapack_freal_state mplapack_freal_tls MPLAPACK_FREAL_TLS_MODEL = {0, 0, nullptr};
std::atomic<unsigned long> mplapack_freal_orphan(0);
std::atomic<unsigned long> mplapack_freal_done(0);

namespace {
long freal_blas_env(const char *name, long def)
{
    const char *s = std::getenv(name);
    return (s && *s) ? std::atol(s) : def;
}
} // namespace

long mplapack_freal_min_work()
{
    static const long w = freal_blas_env("MPLAPACK_MPFR_CPU_FREAL_MIN", 512);
    return w;
}

long mplapack_freal_l12_prec()
{
    static const long p = freal_blas_env("MPLAPACK_MPFR_CPU_FREAL_L12_PREC", 512);
    return p;
}

#define MPLAPACK_FREAL_NS_DECLARE(name, params, args) bool name##_freal params;
#define MPLAPACK_FREAL_PRECISION_NS(P)                       \
    namespace frb##P {                                       \
    MPLAPACK_FREAL_ROUTINES(MPLAPACK_FREAL_NS_DECLARE)       \
    MPLAPACK_FREAL_LAPACK_ROUTINES(MPLAPACK_FREAL_NS_DECLARE) \
    }
MPFREAL_PRECISIONS(MPLAPACK_FREAL_PRECISION_NS)

namespace {

// the precision of the kernels to use, 0: none
int freal_blas_precision()
{
    static const long mode = freal_blas_env("MPLAPACK_MPFR_CPU_FREAL", 1);
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    static const bool cpu = __builtin_cpu_supports("bmi2") && __builtin_cpu_supports("adx");
#else
    static const bool cpu = true;
#endif
    if (mode != 1 || !cpu)
        return 0;
    const mpfr_prec_t p = mpfrxx::default_precision_bits();
    // Fixed-precision kernels construct the literal 1 without MPFR range
    // rounding, including on the one-element STEQR initialization path.
    if (mpfr_get_emin() > 1 || mpfr_get_emax() < 1)
        return 0;
    if (mpfr_get_default_rounding_mode() != MPFR_RNDN || mpfrxx::default_mpc_rounding_mode() != MPC_RNDNN)
        return 0;
    if (mpfrxx::default_mpc_real_precision_bits() != p || mpfrxx::default_mpc_imag_precision_bits() != p)
        return 0;
    switch (p) {
#define MPLAPACK_FREAL_CASE(P) case P:
        MPFREAL_PRECISIONS(MPLAPACK_FREAL_CASE)
#undef MPLAPACK_FREAL_CASE
        return (int)p;
    default:
        return 0;
    }
}

} // namespace

// every precision of MPFREAL_PRECISIONS (Rgemm_freal.h)
#define MPLAPACK_FREAL_CASES(name, args) \
    case 128:                                       \
        return frb128::name##_freal args; \
    case 192:                                       \
        return frb192::name##_freal args; \
    case 256:                                       \
        return frb256::name##_freal args; \
    case 320:                                       \
        return frb320::name##_freal args; \
    case 384:                                       \
        return frb384::name##_freal args; \
    case 448:                                       \
        return frb448::name##_freal args; \
    case 512:                                       \
        return frb512::name##_freal args; \
    case 576:                                       \
        return frb576::name##_freal args; \
    case 640:                                       \
        return frb640::name##_freal args; \
    case 704:                                       \
        return frb704::name##_freal args; \
    case 768:                                       \
        return frb768::name##_freal args; \
    case 832:                                       \
        return frb832::name##_freal args; \
    case 896:                                       \
        return frb896::name##_freal args; \
    case 960:                                       \
        return frb960::name##_freal args; \
    case 1024:                                       \
        return frb1024::name##_freal args; \
    case 1536:                                       \
        return frb1536::name##_freal args; \
    case 2048:                                       \
        return frb2048::name##_freal args;

#define MPLAPACK_FREAL_ENTRY(name, params, args) \
    bool name##_freal params                    \
    {                                           \
        switch (freal_blas_precision()) {       \
            MPLAPACK_FREAL_CASES(name, args)    \
        default:                                \
            return false;                       \
        }                                       \
    }
MPLAPACK_FREAL_ROUTINES(MPLAPACK_FREAL_ENTRY)
MPLAPACK_FREAL_LAPACK_ROUTINES(MPLAPACK_FREAL_ENTRY)
