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

void Ctprfb(const char *side, const char *trans, const char *direct, const char *storev, INTEGER const &m, INTEGER const &n, INTEGER const &k, INTEGER const &l, COMPLEX *v, INTEGER const &ldv, COMPLEX *t, INTEGER const &ldt, COMPLEX *a, INTEGER const &lda, COMPLEX *b, INTEGER const &ldb, COMPLEX *work, INTEGER const &ldwork) {
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
    //     .. Intrinsic Functions ..
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
    const COMPLEX one = (1.0f, 0.0f);
    const COMPLEX zero = (0.0f, 0.0f);
    INTEGER np = 0;
    if (column && forward && left) {
        //
        // ---------------------------------------------------------------------------
        //
        //        Let  W =  [ I ]    (K-by-K)
        //                  [ V ]    (M-by-K)
        //
        //        Form  H C  or  H**H C  where  C = [ A ]  (K-by-N)
        //                                          [ B ]  (M-by-N)
        //
        //        H = I - W T W**H          or  H**H = I - W T**H W**H
        //
        //        A = A -   T (A + V**H B)  or  A = A -   T**H (A + V**H B)
        //        B = B - V T (A + V**H B)  or  B = B - V T**H (A + V**H B)
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
        Ctrmm("L", "U", "C", "N", l, n, one, v[(mp - 1)], ldv, work, ldwork);
        Cgemm("C", "N", l, n, m - l, one, v, ldv, b, ldb, one, work, ldwork);
        Cgemm("C", "N", k - l, n, m, one, v[(kp - 1) * ldv], ldv, b, ldb, zero, work[(kp - 1)], ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("L", "U", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("N", "N", m - l, n, k, -one, v, ldv, work, ldwork, one, b, ldb);
        Cgemm("N", "N", l, n, k - l, -one, v[(mp - 1) + (kp - 1) * ldv], ldv, work[(kp - 1)], ldwork, one, b[(mp - 1)], ldb);
        Ctrmm("L", "U", "N", "N", l, n, one, v[(mp - 1)], ldv, work, ldwork);
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
        //        Form  C H or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N)
        //
        //        H = I - W T W**H          or  H**H = I - W T**H W**H
        //
        //        A = A - (A + B V) T      or  A = A - (A + B V) T**H
        //        B = B - (A + B V) T V**H  or  B = B - (A + B V) T**H V**H
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
        Ctrmm("R", "U", "N", "N", m, l, one, v[(np - 1)], ldv, work, ldwork);
        Cgemm("N", "N", m, l, n - l, one, b, ldb, v, ldv, one, work, ldwork);
        Cgemm("N", "N", m, k - l, n, one, b, ldb, v[(kp - 1) * ldv], ldv, zero, work[(kp - 1) * ldwork], ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("R", "U", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("N", "C", m, n - l, k, -one, work, ldwork, v, ldv, one, b, ldb);
        Cgemm("N", "C", m, l, k - l, -one, work[(kp - 1) * ldwork], ldwork, v[(np - 1) + (kp - 1) * ldv], ldv, one, b[(np - 1) * ldb], ldb);
        Ctrmm("R", "U", "C", "N", m, l, one, v[(np - 1)], ldv, work, ldwork);
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
        //        Form  H C  or  H**H C  where  C = [ B ]  (M-by-N)
        //                                          [ A ]  (K-by-N)
        //
        //        H = I - W T W**H          or  H**H = I - W T**H W**H
        //
        //        A = A -   T (A + V**H B)  or  A = A -   T**H (A + V**H B)
        //        B = B - V T (A + V**H B)  or  B = B - V T**H (A + V**H B)
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
        Ctrmm("L", "L", "C", "N", l, n, one, v[(kp - 1) * ldv], ldv, work[(kp - 1)], ldwork);
        Cgemm("C", "N", l, n, m - l, one, v[(mp - 1) + (kp - 1) * ldv], ldv, b[(mp - 1)], ldb, one, work[(kp - 1)], ldwork);
        Cgemm("C", "N", k - l, n, m, one, v, ldv, b, ldb, zero, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("L", "L", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("N", "N", m - l, n, k, -one, v[(mp - 1)], ldv, work, ldwork, one, b[(mp - 1)], ldb);
        Cgemm("N", "N", l, n, k - l, -one, v, ldv, work, ldwork, one, b, ldb);
        Ctrmm("L", "L", "N", "N", l, n, one, v[(kp - 1) * ldv], ldv, work[(kp - 1)], ldwork);
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
        //        Form  C H  or  C H**H  where  C = [ B A ] (B is M-by-N, A is M-by-K)
        //
        //        H = I - W T W**H          or  H**H = I - W T**H W**H
        //
        //        A = A - (A + B V) T      or  A = A - (A + B V) T**H
        //        B = B - (A + B V) T V**H  or  B = B - (A + B V) T**H V**H
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
        Ctrmm("R", "L", "N", "N", m, l, one, v[(kp - 1) * ldv], ldv, work[(kp - 1) * ldwork], ldwork);
        Cgemm("N", "N", m, l, n - l, one, b[(np - 1) * ldb], ldb, v[(np - 1) + (kp - 1) * ldv], ldv, one, work[(kp - 1) * ldwork], ldwork);
        Cgemm("N", "N", m, k - l, n, one, b, ldb, v, ldv, zero, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("R", "L", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("N", "C", m, n - l, k, -one, work, ldwork, v[(np - 1)], ldv, one, b[(np - 1) * ldb], ldb);
        Cgemm("N", "C", m, l, k - l, -one, work, ldwork, v, ldv, one, b, ldb);
        Ctrmm("R", "L", "C", "N", m, l, one, v[(kp - 1) * ldv], ldv, work[(kp - 1) * ldwork], ldwork);
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
        //        Form  H C  or  H**H C  where  C = [ A ]  (K-by-N)
        //                                          [ B ]  (M-by-N)
        //
        //        H = I - W**H T W          or  H**H = I - W**H T**H W
        //
        //        A = A -     T (A + V B)  or  A = A -     T**H (A + V B)
        //        B = B - V**H T (A + V B)  or  B = B - V**H T**H (A + V B)
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
        Ctrmm("L", "L", "N", "N", l, n, one, v[(mp - 1) * ldv], ldv, work, ldb);
        Cgemm("N", "N", l, n, m - l, one, v, ldv, b, ldb, one, work, ldwork);
        Cgemm("N", "N", k - l, n, m, one, v[(kp - 1)], ldv, b, ldb, zero, work[(kp - 1)], ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("L", "U", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("C", "N", m - l, n, k, -one, v, ldv, work, ldwork, one, b, ldb);
        Cgemm("C", "N", l, n, k - l, -one, v[(kp - 1) + (mp - 1) * ldv], ldv, work[(kp - 1)], ldwork, one, b[(mp - 1)], ldb);
        Ctrmm("L", "L", "C", "N", l, n, one, v[(mp - 1) * ldv], ldv, work, ldwork);
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
        //        Form  C H  or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N)
        //
        //        H = I - W**H T W            or  H**H = I - W**H T**H W
        //
        //        A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H
        //        B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V
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
        Ctrmm("R", "L", "C", "N", m, l, one, v[(np - 1) * ldv], ldv, work, ldwork);
        Cgemm("N", "C", m, l, n - l, one, b, ldb, v, ldv, one, work, ldwork);
        Cgemm("N", "C", m, k - l, n, one, b, ldb, v[(kp - 1)], ldv, zero, work[(kp - 1) * ldwork], ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("R", "U", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("N", "N", m, n - l, k, -one, work, ldwork, v, ldv, one, b, ldb);
        Cgemm("N", "N", m, l, k - l, -one, work[(kp - 1) * ldwork], ldwork, v[(kp - 1) + (np - 1) * ldv], ldv, one, b[(np - 1) * ldb], ldb);
        Ctrmm("R", "L", "N", "N", m, l, one, v[(np - 1) * ldv], ldv, work, ldwork);
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
        //        Form  H C  or  H**H C  where  C = [ B ]  (M-by-N)
        //                                          [ A ]  (K-by-N)
        //
        //        H = I - W**H T W          or  H**H = I - W**H T**H W
        //
        //        A = A -     T (A + V B)  or  A = A -     T**H (A + V B)
        //        B = B - V**H T (A + V B)  or  B = B - V**H T**H (A + V B)
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
        Ctrmm("L", "U", "N", "N", l, n, one, v[(kp - 1)], ldv, work[(kp - 1)], ldwork);
        Cgemm("N", "N", l, n, m - l, one, v[(kp - 1) + (mp - 1) * ldv], ldv, b[(mp - 1)], ldb, one, work[(kp - 1)], ldwork);
        Cgemm("N", "N", k - l, n, m, one, v, ldv, b, ldb, zero, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("L", "L ", trans, "N", k, n, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("C", "N", m - l, n, k, -one, v[(mp - 1) * ldv], ldv, work, ldwork, one, b[(mp - 1)], ldb);
        Cgemm("C", "N", l, n, k - l, -one, v, ldv, work, ldwork, one, b, ldb);
        Ctrmm("L", "U", "C", "N", l, n, one, v[(kp - 1)], ldv, work[(kp - 1)], ldwork);
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
        //        Form  C H  or  C H**H  where  C = [ B A ] (A is M-by-K, B is M-by-N)
        //
        //        H = I - W**H T W            or  H**H = I - W**H T**H W
        //
        //        A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H
        //        B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V
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
        Ctrmm("R", "U", "C", "N", m, l, one, v[(kp - 1)], ldv, work[(kp - 1) * ldwork], ldwork);
        Cgemm("N", "C", m, l, n - l, one, b[(np - 1) * ldb], ldb, v[(kp - 1) + (np - 1) * ldv], ldv, one, work[(kp - 1) * ldwork], ldwork);
        Cgemm("N", "C", m, k - l, n, one, b, ldb, v, ldv, zero, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                work[(i - 1) + (j - 1) * ldwork] += a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        Ctrmm("R", "L", trans, "N", m, k, one, t, ldt, work, ldwork);
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        Cgemm("N", "N", m, n - l, k, -one, work, ldwork, v[(np - 1) * ldv], ldv, one, b[(np - 1) * ldb], ldb);
        Cgemm("N", "N", m, l, k - l, -one, work, ldwork, v, ldv, one, b, ldb);
        Ctrmm("R", "U", "N", "N", m, l, one, v[(kp - 1)], ldv, work[(kp - 1) * ldwork], ldwork);
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - work[(i - 1) + ((k - l + j) - 1) * ldwork];
            }
        }
        //
    }
    //
    //     End of Ctprfb
    //
}
