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

#include <mpblas.h>
#include <mplapack.h>

INTEGER Mmaxloc(REAL const *dx, INTEGER const start, INTEGER const end, INTEGER const incx) {
    // Return the 1-based location of the maximum element within the section:
    //   dx(start:end:incx)
    // Assumptions:
    //   - start/end are 1-based indices into dx
    //   - incx is the stride in that same index space
    // For an empty section or invalid stride (incx == 0), return 0 (Fortran semantics).
    if (dx == nullptr) {
        Mxerbla("Mmaxloc", 1);
        return 0;
    }
    if (incx == 0) {
        Mxerbla("Mmaxloc", 4);
        return 0;
    }

    // Determine whether the section is empty and compute its length n.
    INTEGER n = 0;
    if (incx > 0) {
        if (start > end) {
            return 0;
        }
        n = (end - start) / incx + 1;
    } else { // incx < 0
        if (start < end) {
            return 0;
        }
        n = (start - end) / (-incx) + 1;
    }

    // best_loc: 1-based location in the *section* (Fortran MAXLOC result for rank-1 with DIM=1)
    INTEGER best_loc = 1;
    // best_ix: 1-based index into dx (original index space)
    INTEGER best_ix = start;

    // Fast path for contiguous positive stride.
    if (incx == 1) {
        for (INTEGER loc = 2; loc <= n; ++loc) {
            INTEGER ix = start + (loc - 1);
            // Tie-breaking: keep the first occurrence (Fortran MAXLOC behavior)
            if (dx[ix - 1] > dx[best_ix - 1]) {
                best_ix = ix;
                best_loc = loc;
            }
        }
        return best_loc;
    }

    // General path for arbitrary stride (including negative).
    INTEGER ix = start;
    for (INTEGER loc = 2; loc <= n; ++loc) {
        ix += incx;
        // Tie-breaking: keep the first occurrence
        if (dx[ix - 1] > dx[best_ix - 1]) {
            best_ix = ix;
            best_loc = loc;
        }
    }
    return best_loc;
}
