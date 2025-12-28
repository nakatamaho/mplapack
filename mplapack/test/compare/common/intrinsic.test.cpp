/*
 * Copyright (c) 2008-2025
 *	Nakata, Maho
 * 	All rights reserved.
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
 * Test program for MPLAPACK slice intrinsics helpers:
 *   - Mmaxval
 *   - Mminval
 *   - Mmaxloc
 *
 * This test is modeled after existing *.test.cpp style (e.g., Clarf.test.cpp).
 *
 * Assumed semantics (Fortran-like for rank-1 array sections dx(start:end:incx)):
 *   - start/end are 1-based indices into dx
 *   - incx may be negative (reverse traversal)
 *   - empty section:
 *       (incx > 0 && start > end) || (incx < 0 && start < end)
 *     returns:
 *       Mmaxval -> -Rlamch("O")
 *       Mminval ->  Rlamch("O")
 *       Mmaxloc -> 0
 *   - Mmaxloc returns the position within the SECTION (1..len), not the absolute index.
 *     Ties choose the first occurrence in section order.
 */

#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_compare_debug.h>

#include <cstdio>
#include <cstdlib>
#include <type_traits>

#define MIN_N 1
#define MAX_N 50
#define MAX_ITER 10
#define MAX_ABS_INCX 3
#define MAX_PRINT_ELEMS 12

static inline bool is_empty_section(INTEGER start, INTEGER end, INTEGER incx) {
    if (incx == 0)
        return true;
    if (incx > 0)
        return (start > end);
    return (start < end);
}

static REAL ref_Mmaxval(REAL *dx, INTEGER start, INTEGER end, INTEGER incx) {
    REAL huge = Rlamch("O");
    if (dx == nullptr)
        return -huge;
    if (is_empty_section(start, end, incx))
        return -huge;

    REAL best = dx[start - 1];
    for (INTEGER i = start;; i += incx) {
        if (incx > 0) {
            if (i > end)
                break;
        } else {
            if (i < end)
                break;
        }
        REAL v = dx[i - 1];
        if (v > best)
            best = v;
    }
    return best;
}

static REAL ref_Mminval(REAL *dx, INTEGER start, INTEGER end, INTEGER incx) {
    REAL huge = Rlamch("O");
    if (dx == nullptr)
        return huge;
    if (is_empty_section(start, end, incx))
        return huge;

    REAL best = dx[start - 1];
    for (INTEGER i = start;; i += incx) {
        if (incx > 0) {
            if (i > end)
                break;
        } else {
            if (i < end)
                break;
        }
        REAL v = dx[i - 1];
        if (v < best)
            best = v;
    }
    return best;
}

// Reference maxloc: returns position within the section (1..len), or 0 if empty.
static INTEGER ref_Mmaxloc(REAL *dx, INTEGER start, INTEGER end, INTEGER incx) {
    if (dx == nullptr)
        return 0;
    if (is_empty_section(start, end, incx))
        return 0;

    INTEGER best_pos = 1;
    REAL best = dx[start - 1];

    INTEGER pos = 0;
    for (INTEGER i = start;; i += incx) {
        if (incx > 0) {
            if (i > end)
                break;
        } else {
            if (i < end)
                break;
        }
        pos++;
        REAL v = dx[i - 1];
        // Tie-break: keep the first occurrence (strict > only).
        if (v > best) {
            best = v;
            best_pos = pos;
        }
    }
    return best_pos;
}

static void print_section(REAL *dx, INTEGER n, INTEGER start, INTEGER end, INTEGER incx) {
    (void)n;
    printf("section: [ ");
    int printed = 0;
    for (INTEGER i = start;; i += incx) {
        if (incx > 0) {
            if (i > end)
                break;
        } else {
            if (i < end)
                break;
        }
        if (printed >= MAX_PRINT_ELEMS) {
            printf("... ");
            break;
        }
        printnum(dx[i - 1]);
        printf(" ");
        printed++;
    }
    printf("]\n");
}

// MPLAPACK provides rank-1 maxloc helper:
//   Mmaxloc(dx, start, end, incx)
static INTEGER call_Mmaxloc(REAL *dx, INTEGER start, INTEGER end, INTEGER incx) { return Mmaxloc(dx, start, end, incx); }

static void deterministic_tests() {
    int errorflag = FALSE;

    // Simple vector for deterministic checks.
    // Indices (1-based): 1   2   3   4   5
    // Values:           3   1   3   2   0
    REAL x[5];
    x[0] = (REAL)3;
    x[1] = (REAL)1;
    x[2] = (REAL)3;
    x[3] = (REAL)2;
    x[4] = (REAL)0;

    // Case 1: full range incx=1
    {
        REAL emax = ref_Mmaxval(x, 1, 5, 1);
        REAL emin = ref_Mminval(x, 1, 5, 1);
        INTEGER eloc = ref_Mmaxloc(x, 1, 5, 1);

        REAL gmax = Mmaxval(x, 1, 5, 1);
        REAL gmin = Mminval(x, 1, 5, 1);

        INTEGER gloc = call_Mmaxloc(x, 1, 5, 1);
        if (!(gmax == emax) || !(gmin == emin) || (gloc != eloc)) {
            printf("error: deterministic case 1 failed\n");
            printf("  start=1 end=5 incx=1\n");
            printf("  expected max=");
            printnum(emax);
            printf(" got=");
            printnum(gmax);
            printf("\n");
            printf("  expected min=");
            printnum(emin);
            printf(" got=");
            printnum(gmin);
            printf("\n");
            printf("  expected loc=%d got=%d\n", (int)eloc, (int)gloc);
            errorflag = TRUE;
        }
    }

    // Case 2: tie-breaking (max=3 occurs at positions 1 and 3, expect loc=1)
    {
        INTEGER eloc = ref_Mmaxloc(x, 1, 5, 1);
        INTEGER gloc = call_Mmaxloc(x, 1, 5, 1);
        if (gloc != eloc) {
            printf("error: deterministic tie-break failed\n");
            printf("  expected loc=%d got=%d\n", (int)eloc, (int)gloc);
            errorflag = TRUE;
        }
    }

    // Case 3: strided slice start=1 end=5 incx=2 => elements [3,3,0], max=3, loc=1
    {
        REAL emax = ref_Mmaxval(x, 1, 5, 2);
        REAL emin = ref_Mminval(x, 1, 5, 2);
        INTEGER eloc = ref_Mmaxloc(x, 1, 5, 2);

        REAL gmax = Mmaxval(x, 1, 5, 2);
        REAL gmin = Mminval(x, 1, 5, 2);

        INTEGER gloc = call_Mmaxloc(x, 1, 5, 2);

        if (!(gmax == emax) || !(gmin == emin) || (gloc != eloc)) {
            printf("error: deterministic stride case failed\n");
            printf("  start=1 end=5 incx=2\n");
            printf("  expected max=");
            printnum(emax);
            printf(" got=");
            printnum(gmax);
            printf("\n");
            printf("  expected min=");
            printnum(emin);
            printf(" got=");
            printnum(gmin);
            printf("\n");
            printf("  expected loc=%d got=%d\n", (int)eloc, (int)gloc);
            errorflag = TRUE;
        }
    }

    // Case 4: empty slice (incx>0, start>end)
    {
        REAL huge = Rlamch("O");
        REAL emax = -huge;
        REAL emin = huge;
        INTEGER eloc = 0;

        REAL gmax = Mmaxval(x, 4, 2, 1);
        REAL gmin = Mminval(x, 4, 2, 1);

        INTEGER gloc = call_Mmaxloc(x, 4, 2, 1);
        if (!(gmax == emax) || !(gmin == emin) || (gloc != eloc)) {
            printf("error: deterministic empty case failed\n");
            printf("  start=4 end=2 incx=1\n");
            printf("  expected max=");
            printnum(emax);
            printf(" got=");
            printnum(gmax);
            printf("\n");
            printf("  expected min=");
            printnum(emin);
            printf(" got=");
            printnum(gmin);
            printf("\n");
            printf("  expected loc=%d got=%d\n", (int)eloc, (int)gloc);
            errorflag = TRUE;
        }
    }

    if (errorflag == TRUE) {
        printf("*** Deterministic tests failed ***\n");
        exit(1);
    }
}

static void randomized_tests() {
    int errorflag = FALSE;
    INTEGER n, iter;

    for (n = MIN_N; n <= MAX_N; n++) {
        REAL_REF *x_ref = new REAL_REF[n];
        REAL *x = new REAL[n];

        for (iter = 0; iter < MAX_ITER; iter++) {
            // Fill with random values (same helper style as other tests).
            set_random_vector(x_ref, x, (int)n);

            for (INTEGER incx = -MAX_ABS_INCX; incx <= MAX_ABS_INCX; incx++) {
                if (incx == 0)
                    continue;

                for (INTEGER start = 1; start <= n; start++) {
                    for (INTEGER end = 1; end <= n; end++) {
                        REAL emax = ref_Mmaxval(x, start, end, incx);
                        REAL emin = ref_Mminval(x, start, end, incx);
                        INTEGER eloc = ref_Mmaxloc(x, start, end, incx);

                        REAL gmax = Mmaxval(x, start, end, incx);
                        REAL gmin = Mminval(x, start, end, incx);

                        INTEGER gloc = call_Mmaxloc(x, start, end, incx);
                        bool ok = (gmax == emax) && (gmin == emin) && (gloc == eloc);

                        if (!ok) {
                            printf("error: randomized test failed\n");
                            printf("  n=%d iter=%d start=%d end=%d incx=%d\n", (int)n, (int)iter, (int)start, (int)end, (int)incx);

                            printf("  expected max=");
                            printnum(emax);
                            printf(" got=");
                            printnum(gmax);
                            printf("\n");
                            printf("  expected min=");
                            printnum(emin);
                            printf(" got=");
                            printnum(gmin);
                            printf("\n");
                            printf("  note: Mmaxloc(incx,dim) not available; loc test skipped for incx!=1\n");

                            print_section(x, n, start, end, incx);
                            errorflag = TRUE;
                            break;
                        }
                    }
                    if (errorflag)
                        break;
                }
                if (errorflag)
                    break;
            }
            if (errorflag)
                break;
        }

        delete[] x_ref;
        delete[] x;

        if (errorflag)
            break;
    }

    if (errorflag == TRUE) {
        printf("*** Randomized tests failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    printf("*** Testing Mmaxval/Mminval/Mmaxloc start ***\n");
    deterministic_tests();
    randomized_tests();
    printf("*** Testing Mmaxval/Mminval/Mmaxloc successful ***\n");
    return 0;
}
