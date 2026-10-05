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
 * SIMD kernels for the blocked Rgemm engine (Rgemm_blocked_common.h).
 *
 * dd_simd_kernel: double-double, bit-identical to dd_real (IEEE add).  The
 *   packed planes are the high and low words; the loop over i is vectorized,
 *   each lane performing the dd_real operations of one element of C.
 * qd_bf_kernel: quad-double with Kouya's branch-free addition and
 *   multiplication (qd_real::bf_add, qd_real::bf_mul), vectorized the same
 *   way.  libqd's default qd_real arithmetic (IEEE add) branches on the data
 *   and is not vectorized; this kernel is therefore an opt-in alternative
 *   whose result differs from the default one in the last bits.
 *
 * The block functions are compiled for several instruction sets
 * (Rgemm_simd_common.h); the operations, and hence the results, are the same
 * for every variant.
 */

#ifndef _MPLAPACK_RGEMM_SIMD_KERNELS_COMMON_H_
#define _MPLAPACK_RGEMM_SIMD_KERNELS_COMMON_H_

#include "Rgemm_qd_ops_common.h"
#include "Rgemm_simd_common.h"

namespace mplapack_gemm {

MPLAPACK_GEMM_CLONES static void dd_simd_block(bool axpy, long mc, long nc, long kc, double *acc, long as, const double *a, long ps, const double *b, long bs) {
    for (long j = 0; j < nc; j++) {
        double *__restrict c0 = acc + j * mc;
        double *__restrict c1 = c0 + as;
        for (long l = 0; l < kc; l++) {
            const double *__restrict a0 = a + l * mc;
            const double *__restrict a1 = a0 + ps;
            const double b0 = b[l + j * kc], b1 = b[bs + l + j * kc];
            if (axpy) {
                MPLAPACK_GEMM_SIMD
                for (long i = 0; i < mc; i++) {
                    double p0, p1;
                    dd_mul_w(b0, b1, a0[i], a1[i], p0, p1);
                    dd_add_w(c0[i], c1[i], p0, p1, c0[i], c1[i]);
                }
            } else {
                MPLAPACK_GEMM_SIMD
                for (long i = 0; i < mc; i++) {
                    double p0, p1;
                    dd_mul_w(a0[i], a1[i], b0, b1, p0, p1);
                    dd_add_w(c0[i], c1[i], p0, p1, c0[i], c1[i]);
                }
            }
        }
    }
}

MPLAPACK_GEMM_CLONES static void qd_bf_block(bool axpy, long mc, long nc, long kc, double *acc, long as, const double *a, long ps, const double *b, long bs) {
    for (long j = 0; j < nc; j++) {
        double *__restrict c0 = acc + j * mc;
        double *__restrict c1 = c0 + as;
        double *__restrict c2 = c1 + as;
        double *__restrict c3 = c2 + as;
        for (long l = 0; l < kc; l++) {
            const double *__restrict a0 = a + l * mc;
            const double *__restrict a1 = a0 + ps;
            const double *__restrict a2 = a1 + ps;
            const double *__restrict a3 = a2 + ps;
            const double b0 = b[l + j * kc], b1 = b[bs + l + j * kc], b2 = b[2 * bs + l + j * kc], b3 = b[3 * bs + l + j * kc];
            if (axpy) {
                MPLAPACK_GEMM_SIMD
                for (long i = 0; i < mc; i++) {
                    double p0, p1, p2, p3;
                    qd_bf_mul_w(b0, b1, b2, b3, a0[i], a1[i], a2[i], a3[i], p0, p1, p2, p3);
                    qd_bf_add_w(c0[i], c1[i], c2[i], c3[i], p0, p1, p2, p3, c0[i], c1[i], c2[i], c3[i]);
                }
            } else {
                MPLAPACK_GEMM_SIMD
                for (long i = 0; i < mc; i++) {
                    double p0, p1, p2, p3;
                    qd_bf_mul_w(a0[i], a1[i], a2[i], a3[i], b0, b1, b2, b3, p0, p1, p2, p3);
                    qd_bf_add_w(c0[i], c1[i], c2[i], c3[i], p0, p1, p2, p3, c0[i], c1[i], c2[i], c3[i]);
                }
            }
        }
    }
}

// T is dd_real (two doubles x[0], x[1]).
template <class T> struct dd_simd_kernel {
    typedef T value_type;
    static const int planes = 2;
    static const long MC = 64, NC = 64, KC = 256;
    static void pack(double *p, long stride, long e, const T &v) {
        p[e] = v.x[0];
        p[stride + e] = v.x[1];
    }
    static T unpack(const double *p, long stride, long e) { return T(p[e], p[stride + e]); }
    static void block(bool axpy, long mc, long nc, long kc, double *acc, long as, const double *a, long ps, const double *b, long bs) { dd_simd_block(axpy, mc, nc, kc, acc, as, a, ps, b, bs); }
};

// T is qd_real (four doubles x[0..3]).
template <class T> struct qd_bf_kernel {
    typedef T value_type;
    static const int planes = 4;
    static const long MC = 64, NC = 32, KC = 128;
    static void pack(double *p, long stride, long e, const T &v) {
        for (int q = 0; q < 4; q++)
            p[q * stride + e] = v.x[q];
    }
    static T unpack(const double *p, long stride, long e) { return T(p[e], p[stride + e], p[2 * stride + e], p[3 * stride + e]); }
    static void block(bool axpy, long mc, long nc, long kc, double *acc, long as, const double *a, long ps, const double *b, long bs) { qd_bf_block(axpy, mc, nc, kc, acc, as, a, ps, b, bs); }
};

} // namespace mplapack_gemm

#endif
