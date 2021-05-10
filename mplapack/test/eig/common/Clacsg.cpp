/*
 * Copyright (c) 2021
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

#include <mpblas.h>
#include <mplapack.h>

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

void Clacsg(INTEGER const &m, INTEGER const &p, INTEGER const &q, REAL *theta, INTEGER *iseed, COMPLEX *x, INTEGER const &ldx, COMPLEX *work) {
    INTEGER r = min({p, m - p, q, m - q});
    //
    const COMPLEX zero = COMPLEX(0.0e0, 0.0e0);
    Claset("Full", m, m, zero, zero, x, ldx);
    //
    INTEGER i = 0;
    const COMPLEX one = COMPLEX(1.0e0, 0.0e0);
    for (i = 1; i <= min(p, q) - r; i = i + 1) {
        x[(i - 1) + (i - 1) * ldx] = one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[min(p, q) - r + i - 1 + (min(p, q) - r + i - 1) * ldx] = COMPLEX(cos(theta[i - 1]), 0.0e0);
    }
    for (i = 1; i <= min(p, m - q) - r; i = i + 1) {
        x[((p - i + 1) - 1) + ((m - i + 1) - 1) * ldx] = -one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[p - (min(p, m - q) - r) + 1 - i - 1 + (m - (min(p, m - q) - r) + 1 - i - 1) * ldx] = COMPLEX(-sin(theta[(r - i + 1) - 1]), 0.0e0);
    }
    for (i = 1; i <= min(m - p, q) - r; i = i + 1) {
        x[((m - i + 1) - 1) + ((q - i + 1) - 1) * ldx] = one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[m - (min(m - p, q) - r) + 1 - i - 1 + (q - (min(m - p, q) - r) + 1 - i - 1) * ldx] = COMPLEX(sin(theta[(r - i + 1) - 1]), 0.0e0);
    }
    for (i = 1; i <= min(m - p, m - q) - r; i = i + 1) {
        x[((p + i) - 1) + ((q + i) - 1) * ldx] = one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[p + (min(m - p, m - q) - r) + i - 1 + (q + (min(m - p, m - q) - r) + i - 1) * ldx] = COMPLEX(cos(theta[i - 1]), 0.0e0);
    }
    INTEGER info = 0;
    Claror("Left", "No init", p, m, x, ldx, iseed, work, info);
    Claror("Left", "No init", m - p, m, &x[((p + 1) - 1)], ldx, iseed, work, info);
    Claror("Right", "No init", m, q, x, ldx, iseed, work, info);
    Claror("Right", "No init", m, m - q, &x[((q + 1) - 1) * ldx], ldx, iseed, work, info);
    //
}
