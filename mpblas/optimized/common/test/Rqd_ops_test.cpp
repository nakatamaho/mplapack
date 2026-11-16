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

// The dd/td/qd arithmetic of ../Rgemm_qd_ops_common.h (used by the SIMD and
// GPU kernels) against the operators of libQD3 built with QD_BF: addition,
// subtraction, multiplication, the mixed double forms and the division,
// bit for bit (NaNs aside), on random operands including cancellation, special values,
// huge numerators (rescaled division) and overflow.

#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <limits>
#include <random>
#include "../Rgemm_qd_ops_common.h"

using namespace mplapack_gemm;

namespace {

std::mt19937_64 rng(20261022);
long failures = 0, checks = 0;

double ud(int e) { return std::ldexp(std::uniform_real_distribution<double>(-1.0, 1.0)(rng), e); }

double special()
{
    switch (rng() % 8) {
    case 0: return 0.0;
    case 1: return -0.0;
    case 2: return std::numeric_limits<double>::infinity();
    case 3: return -std::numeric_limits<double>::infinity();
    case 4: return std::numeric_limits<double>::quiet_NaN();
    case 5: return ud(1000);
    case 6: return ud(-1000);
    default: return ud(1023);
    }
}

// Words are compared bit for bit, except that any two NaNs are equal: the
// sign and payload of a NaN depend on the operand order the compiler picks
// for a commutative operation.
void check(const char *what, const double *ref, const double *got, int n)
{
    checks++;
    for (int q = 0; q < n; q++)
        if (std::memcmp(ref + q, got + q, sizeof(double)) != 0 && !(std::isnan(ref[q]) && std::isnan(got[q]))) {
            if (failures++ < 10)
                std::printf("FAIL %s\n", what);
            return;
        }
}

} // namespace

int main()
{
    const long count = 200000;
    for (long t = 0; t < count; t++) {
        double a[4], b[4];
        int ea = (int)(rng() % 400) - 200, eb = ea + (int)(rng() % 9) - 4;
        if (t % 10 == 0)
            ea = 970 + (int)(rng() % 50);
        if (t % 13 == 0)
            eb = -1000 + (int)(rng() % 30);
        for (int q = 0; q < 4; q++)
            a[q] = ud(ea - 53 * q), b[q] = ud(eb - 53 * q);
        if (t % 5 == 1) // cancellation
            for (int q = 0; q < 4; q++)
                b[q] = -a[q] + (q == 3 ? ud(ea - 170) : 0.0);
        const bool sp = t % 7 == 0;
        if (sp)
            a[0] = special(), b[0] = special();
        const double d = b[1];

        dd_real A = dd_real(a[0]), B = dd_real(b[0]);
        if (!sp)
            A = A + dd_real(a[1]), B = B + dd_real(b[1]);
        dd2 x, y, ny, g;
        std::memcpy(x.x, A.x, sizeof x.x), std::memcpy(y.x, B.x, sizeof y.x);
        ny.x[0] = -y.x[0], ny.x[1] = -y.x[1];
        dd_real r;
        r = A + B, g = dd_bf_add(x, y), check("dd +", r.x, g.x, 2);
        r = A - B, g = dd_bf_add(x, ny), check("dd -", r.x, g.x, 2);
        r = A * B, g = dd_bf_mul(x, y), check("dd *", r.x, g.x, 2);
        r = A + d, g = dd_add_d(x, d), check("dd + double", r.x, g.x, 2);
        r = A * d, g = dd_mul_d(x, d), check("dd * double", r.x, g.x, 2);
        r = A / B, g = dd_div(x, y), check("dd /", r.x, g.x, 2);
        r = dd_real::sloppy_add(A, B), g = dd_sloppy_add(x, y), check("dd sloppy +", r.x, g.x, 2);

        td_real TA(a[0], 0.0, 0.0), TB(b[0], 0.0, 0.0);
        if (!sp)
            TA = TA + td_real(a[1]) + td_real(a[2]), TB = TB + td_real(b[1]) + td_real(b[2]);
        const td3 tx = td_make(TA[0], TA[1], TA[2]), ty = td_make(TB[0], TB[1], TB[2]);
        td_real tr;
        td3 tg;
        double tref[3];
#define TD_CHECK(what, expr, port) tr = expr, tg = port, tref[0] = tr[0], tref[1] = tr[1], tref[2] = tr[2], check(what, tref, tg.x, 3)
        TD_CHECK("td +", TA + TB, td_bf_add(tx, ty));
        TD_CHECK("td -", TA - TB, td_bf_add(tx, td_neg(ty)));
        TD_CHECK("td *", TA * TB, td_bf_mul(tx, ty));
        TD_CHECK("td + double", TA + d, td_add_d(tx, d));
        TD_CHECK("td * double", TA * d, td_mul_d(tx, d));
        TD_CHECK("td /", TA / TB, td_div(tx, ty));
        TD_CHECK("td sloppy +", td_real::sloppy_add(TA, TB), td_sloppy_add(tx, ty));

        qd_real QA(a[0]), QB(b[0]);
        if (!sp)
            for (int q = 1; q < 4; q++)
                QA = QA + qd_real(a[q]), QB = QB + qd_real(b[q]);
        const qd4 qx = qd_make(QA[0], QA[1], QA[2], QA[3]), qy = qd_make(QB[0], QB[1], QB[2], QB[3]);
        qd_real qr;
        qd4 qg;
        double qref[4];
#define QD_CHECK(what, expr, port) qr = expr, qg = port, qref[0] = qr[0], qref[1] = qr[1], qref[2] = qr[2], qref[3] = qr[3], check(what, qref, qg.x, 4)
        QD_CHECK("qd +", QA + QB, qd_bf_add(qx, qy));
        QD_CHECK("qd -", QA - QB, qd_bf_add(qx, qd_neg(qy)));
        QD_CHECK("qd *", QA * QB, qd_bf_mul(qx, qy));
        QD_CHECK("qd + double", QA + d, qd_add_d(qx, d));
        QD_CHECK("qd * double", QA * d, qd_mul_d(qx, d));
        QD_CHECK("qd /", QA / QB, qd_div(qx, qy));
        QD_CHECK("qd sloppy +", qd_real::sloppy_add(QA, QB), qd_sloppy_add(qx, qy));
        QD_CHECK("qd sloppy *", qd_real::sloppy_mul(QA, QB), qd_sloppy_mul(qx, qy));
    }
    std::printf("dd/td/qd BF and sloppy arithmetic: %ld checks, %ld failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
