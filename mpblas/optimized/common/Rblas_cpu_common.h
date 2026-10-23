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
 * CPU execution of the lane BLAS (Rblas_lanes_common.h) with OpenMP, for any
 * MPLAPACK type T: the engine's operations are the operators of T, so the
 * result is the reference result bit for bit; the lanes only run in
 * parallel.  MPFR and GMP worker threads get the caller's default precision
 * (and, for MPFR, rounding mode and exponent range), which are thread-local.
 *
 *   MPLAPACK_BLAS_OMP=0         run the reference code (serial)
 *   MPLAPACK_BLAS_OMP_MIN=N     smallest work (multiply-adds) run in parallel
 *                               (default 32768; Level 1: 4 * N)
 *
 * dd/td/qd run the routines of Rblas_simd_common.h instead where it has
 * them (AVX2, branch-free arithmetic; also on a single thread).
 *
 * blas_cpu() returns false, and the caller then runs the reference code, for
 * small calls, OpenMP disabled or a single thread.
 */

#ifndef _MPLAPACK_RBLAS_CPU_COMMON_H_
#define _MPLAPACK_RBLAS_CPU_COMMON_H_

#include <cstdlib>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Rblas_lanes_common.h"
#include "Rblas_calls_common.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "Rblas_simd_common.h"
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
#include <mplapack_mpfr_precision.h>
#endif
#include "Rblas_block_common.h"

namespace mplapack_blas {

template <class T> struct op_engine {
    typedef T V;
    struct M {
        T *p;
        long ld;
    };
    static const bool has_div = true;
    V reg(int) const { return V(0.0); }
    const V &at(const M &a, long i, long j) const { return a.p[i + j * a.ld]; }
    void put(M &a, long i, long j, const V &v) const { a.p[i + j * a.ld] = v; }
    void set(V &r, const V &a) const { r = a; }
    void zero(V &r) const { r = 0.0; }
    void one(V &r) const { r = 1.0; }
    void add(V &r, const V &a, const V &b) const { r = a + b; }
    void sub(V &r, const V &a, const V &b) const { r = a - b; }
    void mul(V &r, const V &a, const V &b) const { r = a * b; }
    void div(V &r, const V &a, const V &b) const { r = a / b; }
    bool iszero(const V &a) const { return a == 0.0; }
    bool isone(const V &a) const { return a == 1.0; }
};

// The caller's arithmetic context, installed in each worker thread.
struct blas_thread_context {
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    mpfr_prec_t prec;
    mpfr_rnd_t rnd;
    mpfr_exp_t emin, emax;
    blas_thread_context() : prec(mpfrxx::default_precision_bits()), rnd(mpfr_get_default_rounding_mode()), emin(mpfr_get_emin()), emax(mpfr_get_emax()) {}
    struct scope {
        MplapackMpfrPrecisionScope p;
        mpfr_rnd_t old_rnd;
        mpfr_exp_t old_emin, old_emax;
        explicit scope(const blas_thread_context &c) : p(c.prec), old_rnd(mpfr_get_default_rounding_mode()), old_emin(mpfr_get_emin()), old_emax(mpfr_get_emax())
        {
            mpfr_set_default_rounding_mode(c.rnd);
            mpfr_set_emin(c.emin);
            mpfr_set_emax(c.emax);
        }
        ~scope()
        {
            mpfr_set_emin(old_emin);
            mpfr_set_emax(old_emax);
            mpfr_set_default_rounding_mode(old_rnd);
        }
    };
#elif defined(MPLAPACK_BUILD_WITH_GMP)
    mp_bitcnt_t prec;
    blas_thread_context() : prec(gmpxx::default_mpf_precision_bits()) {}
    struct scope {
        gmpxx::default_mpf_precision_guard g;
        explicit scope(const blas_thread_context &c) : g(c.prec) {}
    };
#else
    struct scope {
        explicit scope(const blas_thread_context &) {}
    };
#endif
};

template <class T> typename op_engine<T>::M blas_cpu_mat(const blas_operand<T> &o)
{
    typename op_engine<T>::M m;
    m.p = o.X;
    m.ld = o.ld;
    return m;
}

template <class T> bool blas_cpu(blas_call<T> &c)
{
    const bool level1 = c.s.op >= L1_SCAL;
    const double min_work = (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768) * (level1 ? 4 : 1);
    if (c.work < min_work)
        return false;
#ifdef _OPENMP
    const bool omp = blas_env_long("MPLAPACK_BLAS_OMP", 1) != 0 && omp_get_max_threads() >= 2;
#else
    const bool omp = false;
#endif
#ifdef MPLAPACK_BLAS_SIMD_AVAILABLE
    // dd/td/qd: columns and row blocks with AVX2 loops (also on one thread)
    if (blas_simd_cpu(c, omp))
        return true;
#endif
#ifndef _OPENMP
    (void)omp;
    return false;
#else
    if (!omp)
        return false;
    typedef op_engine<T> E;
    const long lanes = blas_lanes(c.s);
    const blas_thread_context ctx;
    std::vector<T> copy;
    blas_mats<E> m;
    m.al.p = &c.alpha, m.al.ld = 1;
    m.be.p = &c.beta, m.be.ld = 1;
    m.a = blas_cpu_mat(c.a);
    m.b = blas_cpu_mat(c.b);
    m.c = blas_cpu_mat(c.c);
    if (c.b_is_copy_of_c) {
        // Rtrmv: the lanes read the original vector
        copy.reserve((size_t)c.c.cols);
        for (long i = 0; i < c.c.cols; i++)
            copy.push_back(c.c.X[i * c.c.ld]); // copies keep the precision
        m.b.p = copy.data();
        m.b.ld = 1;
    }
    const bool heavy = c.s.op <= L3_TRSM;
#pragma omp parallel
    {
        const blas_thread_context::scope sc(ctx);
        E e;
        blas_mats<E> x = m;
        if (heavy) {
#pragma omp for schedule(dynamic, 1)
            for (long lane = 0; lane < lanes; lane++)
                blas_lane(e, c.s, x, lane);
        } else {
#pragma omp for schedule(static)
            for (long lane = 0; lane < lanes; lane++)
                blas_lane(e, c.s, x, lane);
        }
    }
    return true;
#endif
}

} // namespace mplapack_blas

#endif
