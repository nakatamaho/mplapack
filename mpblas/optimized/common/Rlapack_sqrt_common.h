/* Device square roots: operation order of bundled libQD3 1.6.0
 * src/{dd,td,qd}_real.cpp. No approximate reciprocal/sqrt instructions.
 * Adapted from libQD3 (BSD 3-Clause License).
 * Copyright (c) 2021, highprecision. All rights reserved.
 * Copyright (c) 2003-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from U.S. Dept. of Energy). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef MPLAPACK_LAPACK_SQRT_COMMON_H
#define MPLAPACK_LAPACK_SQRT_COMMON_H
#include "Rgemm_cuda_ops_common.h"
namespace mplapack_lapack_gpu {
using namespace mplapack_gemm;
template <class Ops> MPLAPACK_GEMM_HD typename Ops::F power2(typename Ops::F a, double b) {
    for (unsigned i = 0; i < sizeof(a.x) / sizeof(double); ++i)
        a.x[i] *= b;
    return a;
}
template <class Ops> struct root;
template <> struct root<dd_ops> {
    MPLAPACK_GEMM_HD static dd2 core(const dd2 &a) {
        const double x = 1.0 / ::sqrt(a.x[0]), ax = a.x[0] * x;
        double err;
        const double p = ax * ax;
#ifdef QD_FMS
        err = MPLAPACK_GEMM_FMA(ax, ax, -p);
#else
        double hi, lo;
        split(ax, hi, lo);
        err = ((hi * hi - p) + 2.0 * hi * lo) + lo * lo;
#endif
        dd2 square;
        square.x[0] = p;
        square.x[1] = err;
        const double correction = dd_ops::sub(a, square).x[0] * (x * 0.5);
        dd2 r;
        r.x[0] = two_sum(ax, correction, r.x[1]);
        return r;
    }
};
template <> struct root<td_ops> {
    MPLAPACK_GEMM_HD static td3 core(const td3 &a) {
        td3 x = td_ops::from_d(::sqrt(a.x[0]));
        for (int i = 0; i < 2; ++i)
            x = td_ops::add(x, td_ops::div(td_ops::sub(a, td_ops::mul(x, x)), td_mul_d(x, 2.0)));
        return x;
    }
};
template <> struct root<qd_ops> {
    MPLAPACK_GEMM_HD static qd4 core(const qd4 &a) {
        qd4 x = qd_ops::from_d(::sqrt(a.x[0]));
        for (int i = 0; i < 10; ++i) {
            const qd4 y = qd_ops::mul(qd_ops::from_d(0.5), qd_ops::add(x, qd_ops::div(a, x)));
            const qd4 d = qd_ops::sub(x, y);
            x = y;
            if (::fabs(((d.x[3] + d.x[2]) + d.x[1]) + d.x[0]) < ::fabs(x.x[0]) * 0x1p-209)
                break;
        }
        return x;
    }
};
template <class Ops> MPLAPACK_GEMM_HD typename Ops::F square_root(const typename Ops::F &a) {
    if (Ops::is_nan(a) || a.x[0] < 0.0 || MPLAPACK_GEMM_ISINF(a.x[0])) {
        typename Ops::F r;
        const double v = (Ops::is_nan(a) || a.x[0] < 0.0) ? NAN : INFINITY;
        for (unsigned i = 0; i < sizeof(r.x) / sizeof(double); ++i)
            r.x[i] = v;
        return r;
    }
    if (Ops::is_zero(a))
        return a;
    auto scaled = a;
    int shifts = 0;
    if (scaled.x[0] > 0x1p+969) {
        do {
            scaled = power2<Ops>(scaled, 0.25);
            ++shifts;
        } while (sizeof(scaled) == sizeof(td3) && scaled.x[0] > 0x1p+969);
    }
    auto result = root<Ops>::core(scaled);
    while (shifts--)
        result = power2<Ops>(result, 2.0);
    return result;
}
} // namespace mplapack_lapack_gpu
#endif
