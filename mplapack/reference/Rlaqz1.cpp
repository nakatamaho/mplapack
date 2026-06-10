/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DLAQZ1.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Rlaqz1(REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const sr1, REAL const sr2, REAL const si, REAL const beta1, REAL const beta2, REAL *v) {
    //
    // Arguments
    //
    // Parameters
    //
    // Local scalars
    //
    // External Functions
    //
    REAL safmin = Rlamch("SAFE MINIMUM");
    const REAL one = 1.0;
    REAL safmax = one / safmin;
    //
    // Calculate first shifted vector
    //
    REAL w[2];
    w[1 - 1] = beta1 * a[0] - sr1 * b[0];
    w[2 - 1] = beta1 * a[(2 - 1)] - sr1 * b[(2 - 1)];
    REAL scale1 = sqrt(abs(w[1 - 1])) * sqrt(abs(w[2 - 1]));
    if (scale1 >= safmin && scale1 <= safmax) {
        w[1 - 1] = w[1 - 1] / scale1;
        w[2 - 1] = w[2 - 1] / scale1;
    }
    //
    // Solve linear system
    //
    w[2 - 1] = w[2 - 1] / b[(2 - 1) + (2 - 1) * ldb];
    w[1 - 1] = (w[1 - 1] - b[(2 - 1) * ldb] * w[2 - 1]) / b[0];
    REAL scale2 = sqrt(abs(w[1 - 1])) * sqrt(abs(w[2 - 1]));
    if (scale2 >= safmin && scale2 <= safmax) {
        w[1 - 1] = w[1 - 1] / scale2;
        w[2 - 1] = w[2 - 1] / scale2;
    }
    //
    // Apply second shift
    //
    v[1 - 1] = beta2 * (a[0] * w[1 - 1] + a[(2 - 1) * lda] * w[2 - 1]) - sr2 * (b[0] * w[1 - 1] + b[(2 - 1) * ldb] * w[2 - 1]);
    v[2 - 1] = beta2 * (a[(2 - 1)] * w[1 - 1] + a[(2 - 1) + (2 - 1) * lda] * w[2 - 1]) - sr2 * (b[(2 - 1)] * w[1 - 1] + b[(2 - 1) + (2 - 1) * ldb] * w[2 - 1]);
    v[3 - 1] = beta2 * (a[(3 - 1)] * w[1 - 1] + a[(3 - 1) + (2 - 1) * lda] * w[2 - 1]) - sr2 * (b[(3 - 1)] * w[1 - 1] + b[(3 - 1) + (2 - 1) * ldb] * w[2 - 1]);
    //
    // Account for imaginary part
    //
    v[1 - 1] += si * si * b[0] / scale1 / scale2;
    //
    // Check for overflow
    //
    const REAL zero = 0.0;
    if (abs(v[1 - 1]) > safmax || abs(v[2 - 1]) > safmax || abs(v[3 - 1]) > safmax || Risnan(v[1 - 1]) || Risnan(v[2 - 1]) || Risnan(v[3 - 1])) {
        v[1 - 1] = zero;
        v[2 - 1] = zero;
        v[3 - 1] = zero;
    }
    //
    // End of Rlaqz1
    //
}
