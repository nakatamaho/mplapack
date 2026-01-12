/*
 * Copyright (c) 2021-2025
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

REAL Mmaxval(REAL *dx, INTEGER const start, INTEGER const end, INTEGER const incx) {
    // Fortran MAXVAL for a rank-1 section:
    //   maxval( dx(start:end:incx) )
    // Assumptions:
    //   - start/end are 1-based indices into dx
    //   - incx is the section stride (may be negative)
    // Semantics for an empty section:
    // - return -HUGE (implemented as -Rlamch("O") in LAPACK/MPLAPACK style)

    if (dx == nullptr) {
        Mxerbla("Mmaxval", 1);
        return -Rlamch("O");
    }

    if (incx == 0) {
        // Invalid stride: treat as empty section (Fortran array section with stride 0 is illegal)
        Mxerbla("Mmaxval", 4);
        return -Rlamch("O");
    }

    // Empty section check
    if ((incx > 0 && start > end) || (incx < 0 && start < end)) {
        return -Rlamch("O");
    }

    REAL dmax = dx[start - 1];

    // Fast path: contiguous forward
    if (incx == 1) {
        for (INTEGER i = start + 1; i <= end; ++i) {
            if (dx[i - 1] > dmax) {
                dmax = dx[i - 1];
            }
        }
        return dmax;
    }

    // General path: arbitrary stride (positive or negative)
    INTEGER ix = start + incx;
    if (incx > 0) {
        for (; ix <= end; ix += incx) {
            if (dx[ix - 1] > dmax) {
                dmax = dx[ix - 1];
            }
        }
    } else { // incx < 0
        for (; ix >= end; ix += incx) {
            if (dx[ix - 1] > dmax) {
                dmax = dx[ix - 1];
            }
        }
    }

    return dmax;
}
