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

// The per-thread state of the fixed-precision MPFR BLAS (Rfreal_blas.h):
// the exponent range a result must stay in to be the MPFR result, and where
// to report one that leaves it (or an operation the kernels do not have).
// A call that sees a report gives up and the MPFR code runs instead.
// blas_thread_context (../../common/Rblas_cpu_common.h) carries it into the
// OpenMP worker threads; a worker without it has the empty range, so every
// nonzero result there is reported (to mplapack_freal_orphan).

#ifndef _MPLAPACK_RFREAL_STATE_H_
#define _MPLAPACK_RFREAL_STATE_H_

#include <atomic>

struct mplapack_freal_state {
    long emin, emax;
    std::atomic<bool> *bad;
};

#if defined(__GNUC__) && !defined(_WIN32)
#define MPLAPACK_FREAL_TLS_MODEL __attribute__((tls_model("initial-exec")))
#else
#define MPLAPACK_FREAL_TLS_MODEL
#endif

extern thread_local mplapack_freal_state mplapack_freal_tls MPLAPACK_FREAL_TLS_MODEL;
extern std::atomic<unsigned long> mplapack_freal_orphan;
// calls done by the kernels (for the tests)
extern std::atomic<unsigned long> mplapack_freal_done;
// MPLAPACK_MPFR_CPU_FREAL_MIN, read at the first call
long mplapack_freal_min_work();
// MPLAPACK_MPFR_CPU_FREAL_L12_PREC: the largest precision for Level 1 and 2
long mplapack_freal_l12_prec();

inline void mplapack_freal_fail()
{
    if (mplapack_freal_tls.bad)
        mplapack_freal_tls.bad->store(true, std::memory_order_relaxed);
    else
        mplapack_freal_orphan.fetch_add(1, std::memory_order_relaxed);
}

#endif
