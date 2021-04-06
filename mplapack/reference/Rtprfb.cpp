/*
 * Copyright (c) 2008-2021
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

void Rtprfb(const char *side, const char *trans, const char *direct, const char *storev, INTEGER const m, INTEGER const n, INTEGER const k, INTEGER const l, REAL *v, INTEGER const ldv, REAL *t, INTEGER const ldt, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *work, INTEGER const ldwork) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  ==========================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0 || k <= 0 || l < 0) {
        return;
    }
    //
    bool column = false;
    bool row = false;
    if (Mlsame(storev, "C")) {
        column = true;
        row = false;
    } else if (Mlsame(storev, "R")) {
        column = false;
        row = true;
    } else {
        column = false;
        row = false;
    }
    //
    bool left = false;
    bool right = false;
    if (Mlsame(side, "L")) {
        left = true;
        right = false;
    } else if (Mlsame(side, "R")) {
        left = false;
        right = true;
    } else {
        left = false;
        right = false;
    }
    //
    bool forward = false;
    bool backward = false;
    if (Mlsame(direct, "F")) {
        forward = true;
        backward = false;
    } else if (Mlsame(direct, "B")) {
        forward = false;
        backward = true;
    } else {
        forward = false;
        backward = false;
    }
    //
    // ---------------------------------------------------------------------------
    //
    INTEGER mp = 0;
    INTEGER kp = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL one = 1.0f;
    const REAL zero = 0.0f;
    INTEGER np = 0;
    if (column && forward && left) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ I ]    (K-by-K)
        //                  [ V ]    (M-by-K)
        //
        //        Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
        //                                          [ B ]  (M-by-N)
        //
        //        H = I - W T W**T          or  H**T = I - W T**T W**T
        //
        //        A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
        //        B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)
        //
        // ---------------------------------------------------------------------------
        //
        mp = min(m - l + 1, m);
        kp = min(l + 1, k);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] = b[((m - l + i) - 1) + (j - 1) * ldb];
            }
        }
        Rtrmm("L", "U", "T", "N", l, n, one, &v[(mp - 1)], ldv, work, ldwork);
        Rgemm("T", "N", l, n, m - l, one, v, ldv, b, ldb, one, work, ldwork);
        Rgemm("T", "N", k - l, n, m, one, &v[(kp - 1) * ldv], ldv, b, ldb, zero, &work[(kp - 1)], ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("L", "U", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("N", "N", m - l, n, k, -one, v, ldv, work, ldwork, one, b, ldb);
        Rgemm("N", "N", l, n, k - l, -one, &v[(mp - 1) + (kp - 1) * ldv], ldv, &work[(kp - 1)], ldwork, one, &b[(mp - 1)], ldb);
        Rtrmm("L", "U", "N", "N", l, n, one, &v[(mp - 1)], ldv, work, ldwork);
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                b[((m - l + i) - 1) + (j - 1) * ldb] = b[((m - l + i) - 1) + (j - 1) * ldb] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (column && forward && right) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ I ]    (K-by-K)
        //                  [ V ]    (N-by-K)
        //
        //        Form  C H or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)
        //
        //        H = I - W T W**T          or  H**T = I - W T**T W**T
        //
        //        A = A - (A + B V) T      or  A = A - (A + B V) T**T
        //        B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T
        //
        // ---------------------------------------------------------------------------
        //
        np = min(n - l + 1, n);
        kp = min(l + 1, k);
        //
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] = b[(i - 1) + ((n - l + j) - 1) * ldb];
            }
        }
        Rtrmm("R", "U", "N", "N", m, l, one, &v[(np - 1)], ldv, work, ldwork);
        Rgemm("N", "N", m, l, n - l, one, b, ldb, v, ldv, one, work, ldwork);
        Rgemm("N", "N", m, k - l, n, one, b, ldb, &v[(kp - 1) * ldv], ldv, zero, &work[(kp - 1) * ldwork], ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("R", "U", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("N", "T", m, n - l, k, -one, work, ldwork, v, ldv, one, b, ldb);
        Rgemm("N", "T", m, l, k - l, -one, &work[(kp - 1) * ldwork], ldwork, &v[(np - 1) + (kp - 1) * ldv], ldv, one, &b[(np - 1) * ldb], ldb);
        Rtrmm("R", "U", "T", "N", m, l, one, &v[(np - 1)], ldv, work, ldwork);
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + ((n - l + j) - 1) * ldb] = b[(i - 1) + ((n - l + j) - 1) * ldb] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (column && backward && left) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ V ]    (M-by-K)
        //                  [ I ]    (K-by-K)
        //
        //        Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
        //                                          [ A ]  (K-by-N)
        //
        //        H = I - W T W**T          or  H**T = I - W T**T W**T
        //
        //        A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
        //        B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)
        //
        // ---------------------------------------------------------------------------
        //
        mp = min(l + 1, m);
        kp = min(k - l + 1, k);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                work[((k - l + i) - 1) + (j - 1) * ldwork] = b[(i - 1) + (j - 1) * ldb];
            }
        }
        //
        Rtrmm("L", "L", "T", "N", l, n, one, &v[(kp - 1) * ldv], ldv, &work[(kp - 1)], ldwork);
        Rgemm("T", "N", l, n, m - l, one, &v[(mp - 1) + (kp - 1) * ldv], ldv, &b[(mp - 1)], ldb, one, &work[(kp - 1)], ldwork);
        Rgemm("T", "N", k - l, n, m, one, v, ldv, b, ldb, zero, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("L", "L", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("N", "N", m - l, n, k, -one, &v[(mp - 1)], ldv, work, ldwork, one, &b[(mp - 1)], ldb);
        Rgemm("N", "N", l, n, k - l, -one, v, ldv, work, ldwork, one, b, ldb);
        Rtrmm("L", "L", "N", "N", l, n, one, &v[(kp - 1) * ldv], ldv, &work[(kp - 1)], ldwork);
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - work[((k - l + i) - 1) + (j - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (column && backward && right) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ V ]    (N-by-K)
        //                  [ I ]    (K-by-K)
        //
        //        Form  C H  or  C H**T  where  C = [ B A ] (B is M-by-N, A is M-by-K)
        //
        //        H = I - W T W**T          or  H**T = I - W T**T W**T
        //
        //        A = A - (A + B V) T      or  A = A - (A + B V) T**T
        //        B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T
        //
        // ---------------------------------------------------------------------------
        //
        np = min(l + 1, n);
        kp = min(k - l + 1, k);
        //
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + ((k - l + j) - 1) * ldwork] = b[(i - 1) + (j - 1) * ldb];
            }
        }
        Rtrmm("R", "L", "N", "N", m, l, one, &v[(kp - 1) * ldv], ldv, &work[(kp - 1) * ldwork], ldwork);
        Rgemm("N", "N", m, l, n - l, one, &b[(np - 1) * ldb], ldb, &v[(np - 1) + (kp - 1) * ldv], ldv, one, &work[(kp - 1) * ldwork], ldwork);
        Rgemm("N", "N", m, k - l, n, one, b, ldb, v, ldv, zero, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("R", "L", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("N", "T", m, n - l, k, -one, work, ldwork, &v[(np - 1)], ldv, one, &b[(np - 1) * ldb], ldb);
        Rgemm("N", "T", m, l, k - l, -one, work, ldwork, v, ldv, one, b, ldb);
        Rtrmm("R", "L", "T", "N", m, l, one, &v[(kp - 1) * ldv], ldv, &work[(kp - 1) * ldwork], ldwork);
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - work[(i - 1) + ((k - l + j) - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (row && forward && left) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ I V ] ( I is K-by-K, V is K-by-M )
        //
        //        Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
        //                                          [ B ]  (M-by-N)
        //
        //        H = I - W**T T W          or  H**T = I - W**T T**T W
        //
        //        A = A -     T (A + V B)  or  A = A -     T**T (A + V B)
        //        B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)
        //
        // ---------------------------------------------------------------------------
        //
        mp = min(m - l + 1, m);
        kp = min(l + 1, k);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] = b[((m - l + i) - 1) + (j - 1) * ldb];
            }
        }
        Rtrmm("L", "L", "N", "N", l, n, one, &v[(mp - 1) * ldv], ldv, work, ldb);
        Rgemm("N", "N", l, n, m - l, one, v, ldv, b, ldb, one, work, ldwork);
        Rgemm("N", "N", k - l, n, m, one, &v[(kp - 1)], ldv, b, ldb, zero, &work[(kp - 1)], ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("L", "U", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("T", "N", m - l, n, k, -one, v, ldv, work, ldwork, one, b, ldb);
        Rgemm("T", "N", l, n, k - l, -one, &v[(kp - 1) + (mp - 1) * ldv], ldv, &work[(kp - 1)], ldwork, one, &b[(mp - 1)], ldb);
        Rtrmm("L", "L", "T", "N", l, n, one, &v[(mp - 1) * ldv], ldv, work, ldwork);
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                b[((m - l + i) - 1) + (j - 1) * ldb] = b[((m - l + i) - 1) + (j - 1) * ldb] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (row && forward && right) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ I V ] ( I is K-by-K, V is K-by-N )
        //
        //        Form  C H  or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)
        //
        //        H = I - W**T T W            or  H**T = I - W**T T**T W
        //
        //        A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
        //        B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V
        //
        // ---------------------------------------------------------------------------
        //
        np = min(n - l + 1, n);
        kp = min(l + 1, k);
        //
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] = b[(i - 1) + ((n - l + j) - 1) * ldb];
            }
        }
        Rtrmm("R", "L", "T", "N", m, l, one, &v[(np - 1) * ldv], ldv, work, ldwork);
        Rgemm("N", "T", m, l, n - l, one, b, ldb, v, ldv, one, work, ldwork);
        Rgemm("N", "T", m, k - l, n, one, b, ldb, &v[(kp - 1)], ldv, zero, &work[(kp - 1) * ldwork], ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("R", "U", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("N", "N", m, n - l, k, -one, work, ldwork, v, ldv, one, b, ldb);
        Rgemm("N", "N", m, l, k - l, -one, &work[(kp - 1) * ldwork], ldwork, &v[(kp - 1) + (np - 1) * ldv], ldv, one, &b[(np - 1) * ldb], ldb);
        Rtrmm("R", "L", "N", "N", m, l, one, &v[(np - 1) * ldv], ldv, work, ldwork);
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + ((n - l + j) - 1) * ldb] = b[(i - 1) + ((n - l + j) - 1) * ldb] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (row && backward && left) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ V I ] ( I is K-by-K, V is K-by-M )
        //
        //        Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
        //                                          [ A ]  (K-by-N)
        //
        //        H = I - W**T T W          or  H**T = I - W**T T**T W
        //
        //        A = A -     T (A + V B)  or  A = A -     T**T (A + V B)
        //        B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)
        //
        // ---------------------------------------------------------------------------
        //
        mp = min(l + 1, m);
        kp = min(k - l + 1, k);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                work[((k - l + i) - 1) + (j - 1) * ldwork] = b[(i - 1) + (j - 1) * ldb];
            }
        }
        Rtrmm("L", "U", "N", "N", l, n, one, &v[(kp - 1)], ldv, &work[(kp - 1)], ldwork);
        Rgemm("N", "N", l, n, m - l, one, &v[(kp - 1) + (mp - 1) * ldv], ldv, &b[(mp - 1)], ldb, one, &work[(kp - 1)], ldwork);
        Rgemm("N", "N", k - l, n, m, one, v, ldv, b, ldb, zero, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("L", "L ", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("T", "N", m - l, n, k, -one, &v[(mp - 1) * ldv], ldv, work, ldwork, one, &b[(mp - 1)], ldb);
        Rgemm("T", "N", l, n, k - l, -one, v, ldv, work, ldwork, one, b, ldb);
        Rtrmm("L", "U", "T", "N", l, n, one, &v[(kp - 1)], ldv, &work[(kp - 1)], ldwork);
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - work[((k - l + i) - 1) + (j - 1) * ldwork];
            }
        }
        //
        // ---------------------------------------------------------------------------
        //
    } else if (row && backward && right) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ V I ] ( I is K-by-K, V is K-by-N )
        //
        //        Form  C H  or  C H**T  where  C = [ B A ] (A is M-by-K, B is M-by-N)
        //
        //        H = I - W**T T W            or  H**T = I - W**T T**T W
        //
        //        A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
        //        B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V
        //
        // ---------------------------------------------------------------------------
        //
        np = min(l + 1, n);
        kp = min(k - l + 1, k);
        //
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + ((k - l + j) - 1) * ldwork] = b[(i - 1) + (j - 1) * ldb];
            }
        }
        Rtrmm("R", "U", "T", "N", m, l, one, &v[(kp - 1)], ldv, &work[(kp - 1) * ldwork], ldwork);
        Rgemm("N", "T", m, l, n - l, one, &b[(np - 1) * ldb], ldb, &v[(kp - 1) + (np - 1) * ldv], ldv, one, &work[(kp - 1) * ldwork], ldwork);
        Rgemm("N", "T", m, k - l, n, one, b, ldb, v, ldv, zero, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Rtrmm("R", "L", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Rgemm("N", "N", m, n - l, k, -one, work, ldwork, &v[(np - 1) * ldv], ldv, one, &b[(np - 1) * ldb], ldb);
        Rgemm("N", "N", m, l, k - l, -one, work, ldwork, v, ldv, one, b, ldb);
        Rtrmm("R", "U", "N", "N", m, l, one, &v[(kp - 1)], ldv, &work[(kp - 1) * ldwork], ldwork);
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - work[(i - 1) + ((k - l + j) - 1) * ldwork];
            }
        }
        //
    }
    //
    //     End of Rtprfb
    //
}
