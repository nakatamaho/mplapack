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

void Rlarfb(const char *side, const char *trans, const char *direct, const char *storev, INTEGER const &m, INTEGER const &n, INTEGER const &k, REAL *v, INTEGER const &ldv, REAL *t, INTEGER const &ldt, REAL *c, INTEGER const &ldc, REAL *work, INTEGER const &ldwork) {
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
    //  =====================================================================
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
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    str<1> transt = char0;
    if (Mlsame(trans, "N")) {
        transt = "T";
    } else {
        transt = "N";
    }
    //
    INTEGER j = 0;
    const REAL one = 1.0;
    INTEGER i = 0;
    if (Mlsame(storev, "C")) {
        //
        if (Mlsame(direct, "F")) {
            //
            //           Let  V =  ( V1 )    (first K rows)
            //                     ( V2 )
            //           where  V1  is unit lower triangular.
            //
            if (Mlsame(side, "L")) {
                //
                //              Form  H * C  or  H**T * C  where  C = ( C1 )
                //                                                    ( C2 )
                //
                //              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
                //
                //              W := C1**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(n, c[(j - 1)], ldc, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V1
                //
                Rtrmm("Right", "Lower", "No transpose", "Unit", n, k, one, v, ldv, work, ldwork);
                if (m > k) {
                    //
                    //                 W := W + C2**T * V2
                    //
                    Rgemm("Transpose", "No transpose", n, k, m - k, one, c[((k + 1) - 1)], ldc, v[((k + 1) - 1)], ldv, one, work, ldwork);
                }
                //
                //              W := W * T**T  or  W * T
                //
                Rtrmm("Right", "Upper", transt, "Non-unit", n, k, one, t, ldt, work, ldwork);
                //
                //              C := C - V * W**T
                //
                if (m > k) {
                    //
                    //                 C2 := C2 - V2 * W**T
                    //
                    Rgemm("No transpose", "Transpose", m - k, n, k, -one, v[((k + 1) - 1)], ldv, work, ldwork, one, c[((k + 1) - 1)], ldc);
                }
                //
                //              W := W * V1**T
                //
                Rtrmm("Right", "Lower", "Transpose", "Unit", n, k, one, v, ldv, work, ldwork);
                //
                //              C1 := C1 - W**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        c[(j - 1) + (i - 1) * ldc] = c[(j - 1) + (i - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            } else if (Mlsame(side, "R")) {
                //
                //              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
                //
                //              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
                //
                //              W := C1
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(m, c[(j - 1) * ldc], 1, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V1
                //
                Rtrmm("Right", "Lower", "No transpose", "Unit", m, k, one, v, ldv, work, ldwork);
                if (n > k) {
                    //
                    //                 W := W + C2 * V2
                    //
                    Rgemm("No transpose", "No transpose", m, k, n - k, one, c[((k + 1) - 1) * ldc], ldc, v[((k + 1) - 1)], ldv, one, work, ldwork);
                }
                //
                //              W := W * T  or  W * T**T
                //
                Rtrmm("Right", "Upper", trans, "Non-unit", m, k, one, t, ldt, work, ldwork);
                //
                //              C := C - W * V**T
                //
                if (n > k) {
                    //
                    //                 C2 := C2 - W * V2**T
                    //
                    Rgemm("No transpose", "Transpose", m, n - k, k, -one, work, ldwork, v[((k + 1) - 1)], ldv, one, c[((k + 1) - 1) * ldc], ldc);
                }
                //
                //              W := W * V1**T
                //
                Rtrmm("Right", "Lower", "Transpose", "Unit", m, k, one, v, ldv, work, ldwork);
                //
                //              C1 := C1 - W
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
            }
            //
        } else {
            //
            //           Let  V =  ( V1 )
            //                     ( V2 )    (last K rows)
            //           where  V2  is unit upper triangular.
            //
            if (Mlsame(side, "L")) {
                //
                //              Form  H * C  or  H**T * C  where  C = ( C1 )
                //                                                    ( C2 )
                //
                //              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
                //
                //              W := C2**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(n, c[((m - k + j) - 1)], ldc, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V2
                //
                Rtrmm("Right", "Upper", "No transpose", "Unit", n, k, one, v[((m - k + 1) - 1)], ldv, work, ldwork);
                if (m > k) {
                    //
                    //                 W := W + C1**T * V1
                    //
                    Rgemm("Transpose", "No transpose", n, k, m - k, one, c, ldc, v, ldv, one, work, ldwork);
                }
                //
                //              W := W * T**T  or  W * T
                //
                Rtrmm("Right", "Lower", transt, "Non-unit", n, k, one, t, ldt, work, ldwork);
                //
                //              C := C - V * W**T
                //
                if (m > k) {
                    //
                    //                 C1 := C1 - V1 * W**T
                    //
                    Rgemm("No transpose", "Transpose", m - k, n, k, -one, v, ldv, work, ldwork, one, c, ldc);
                }
                //
                //              W := W * V2**T
                //
                Rtrmm("Right", "Upper", "Transpose", "Unit", n, k, one, v[((m - k + 1) - 1)], ldv, work, ldwork);
                //
                //              C2 := C2 - W**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        c[((m - k + j) - 1) + (i - 1) * ldc] = c[((m - k + j) - 1) + (i - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            } else if (Mlsame(side, "R")) {
                //
                //              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
                //
                //              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
                //
                //              W := C2
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(m, c[((n - k + j) - 1) * ldc], 1, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V2
                //
                Rtrmm("Right", "Upper", "No transpose", "Unit", m, k, one, v[((n - k + 1) - 1)], ldv, work, ldwork);
                if (n > k) {
                    //
                    //                 W := W + C1 * V1
                    //
                    Rgemm("No transpose", "No transpose", m, k, n - k, one, c, ldc, v, ldv, one, work, ldwork);
                }
                //
                //              W := W * T  or  W * T**T
                //
                Rtrmm("Right", "Lower", trans, "Non-unit", m, k, one, t, ldt, work, ldwork);
                //
                //              C := C - W * V**T
                //
                if (n > k) {
                    //
                    //                 C1 := C1 - W * V1**T
                    //
                    Rgemm("No transpose", "Transpose", m, n - k, k, -one, work, ldwork, v, ldv, one, c, ldc);
                }
                //
                //              W := W * V2**T
                //
                Rtrmm("Right", "Upper", "Transpose", "Unit", m, k, one, v[((n - k + 1) - 1)], ldv, work, ldwork);
                //
                //              C2 := C2 - W
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + ((n - k + j) - 1) * ldc] = c[(i - 1) + ((n - k + j) - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
            }
        }
        //
    } else if (Mlsame(storev, "R")) {
        //
        if (Mlsame(direct, "F")) {
            //
            //           Let  V =  ( V1  V2 )    (V1: first K columns)
            //           where  V1  is unit upper triangular.
            //
            if (Mlsame(side, "L")) {
                //
                //              Form  H * C  or  H**T * C  where  C = ( C1 )
                //                                                    ( C2 )
                //
                //              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
                //
                //              W := C1**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(n, c[(j - 1)], ldc, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V1**T
                //
                Rtrmm("Right", "Upper", "Transpose", "Unit", n, k, one, v, ldv, work, ldwork);
                if (m > k) {
                    //
                    //                 W := W + C2**T * V2**T
                    //
                    Rgemm("Transpose", "Transpose", n, k, m - k, one, c[((k + 1) - 1)], ldc, v[((k + 1) - 1) * ldv], ldv, one, work, ldwork);
                }
                //
                //              W := W * T**T  or  W * T
                //
                Rtrmm("Right", "Upper", transt, "Non-unit", n, k, one, t, ldt, work, ldwork);
                //
                //              C := C - V**T * W**T
                //
                if (m > k) {
                    //
                    //                 C2 := C2 - V2**T * W**T
                    //
                    Rgemm("Transpose", "Transpose", m - k, n, k, -one, v[((k + 1) - 1) * ldv], ldv, work, ldwork, one, c[((k + 1) - 1)], ldc);
                }
                //
                //              W := W * V1
                //
                Rtrmm("Right", "Upper", "No transpose", "Unit", n, k, one, v, ldv, work, ldwork);
                //
                //              C1 := C1 - W**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        c[(j - 1) + (i - 1) * ldc] = c[(j - 1) + (i - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            } else if (Mlsame(side, "R")) {
                //
                //              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
                //
                //              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
                //
                //              W := C1
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(m, c[(j - 1) * ldc], 1, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V1**T
                //
                Rtrmm("Right", "Upper", "Transpose", "Unit", m, k, one, v, ldv, work, ldwork);
                if (n > k) {
                    //
                    //                 W := W + C2 * V2**T
                    //
                    Rgemm("No transpose", "Transpose", m, k, n - k, one, c[((k + 1) - 1) * ldc], ldc, v[((k + 1) - 1) * ldv], ldv, one, work, ldwork);
                }
                //
                //              W := W * T  or  W * T**T
                //
                Rtrmm("Right", "Upper", trans, "Non-unit", m, k, one, t, ldt, work, ldwork);
                //
                //              C := C - W * V
                //
                if (n > k) {
                    //
                    //                 C2 := C2 - W * V2
                    //
                    Rgemm("No transpose", "No transpose", m, n - k, k, -one, work, ldwork, v[((k + 1) - 1) * ldv], ldv, one, c[((k + 1) - 1) * ldc], ldc);
                }
                //
                //              W := W * V1
                //
                Rtrmm("Right", "Upper", "No transpose", "Unit", m, k, one, v, ldv, work, ldwork);
                //
                //              C1 := C1 - W
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            }
            //
        } else {
            //
            //           Let  V =  ( V1  V2 )    (V2: last K columns)
            //           where  V2  is unit lower triangular.
            //
            if (Mlsame(side, "L")) {
                //
                //              Form  H * C  or  H**T * C  where  C = ( C1 )
                //                                                    ( C2 )
                //
                //              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
                //
                //              W := C2**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(n, c[((m - k + j) - 1)], ldc, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V2**T
                //
                Rtrmm("Right", "Lower", "Transpose", "Unit", n, k, one, v[((m - k + 1) - 1) * ldv], ldv, work, ldwork);
                if (m > k) {
                    //
                    //                 W := W + C1**T * V1**T
                    //
                    Rgemm("Transpose", "Transpose", n, k, m - k, one, c, ldc, v, ldv, one, work, ldwork);
                }
                //
                //              W := W * T**T  or  W * T
                //
                Rtrmm("Right", "Lower", transt, "Non-unit", n, k, one, t, ldt, work, ldwork);
                //
                //              C := C - V**T * W**T
                //
                if (m > k) {
                    //
                    //                 C1 := C1 - V1**T * W**T
                    //
                    Rgemm("Transpose", "Transpose", m - k, n, k, -one, v, ldv, work, ldwork, one, c, ldc);
                }
                //
                //              W := W * V2
                //
                Rtrmm("Right", "Lower", "No transpose", "Unit", n, k, one, v[((m - k + 1) - 1) * ldv], ldv, work, ldwork);
                //
                //              C2 := C2 - W**T
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        c[((m - k + j) - 1) + (i - 1) * ldc] = c[((m - k + j) - 1) + (i - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            } else if (Mlsame(side, "R")) {
                //
                //              Form  C * H  or  C * H'  where  C = ( C1  C2 )
                //
                //              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
                //
                //              W := C2
                //
                for (j = 1; j <= k; j = j + 1) {
                    Rcopy(m, c[((n - k + j) - 1) * ldc], 1, work[(j - 1) * ldwork], 1);
                }
                //
                //              W := W * V2**T
                //
                Rtrmm("Right", "Lower", "Transpose", "Unit", m, k, one, v[((n - k + 1) - 1) * ldv], ldv, work, ldwork);
                if (n > k) {
                    //
                    //                 W := W + C1 * V1**T
                    //
                    Rgemm("No transpose", "Transpose", m, k, n - k, one, c, ldc, v, ldv, one, work, ldwork);
                }
                //
                //              W := W * T  or  W * T**T
                //
                Rtrmm("Right", "Lower", trans, "Non-unit", m, k, one, t, ldt, work, ldwork);
                //
                //              C := C - W * V
                //
                if (n > k) {
                    //
                    //                 C1 := C1 - W * V1
                    //
                    Rgemm("No transpose", "No transpose", m, n - k, k, -one, work, ldwork, v, ldv, one, c, ldc);
                }
                //
                //              W := W * V2
                //
                Rtrmm("Right", "Lower", "No transpose", "Unit", m, k, one, v[((n - k + 1) - 1) * ldv], ldv, work, ldwork);
                //
                //              C1 := C1 - W
                //
                for (j = 1; j <= k; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + ((n - k + j) - 1) * ldc] = c[(i - 1) + ((n - k + j) - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            }
            //
        }
    }
    //
    //     End of Rlarfb
    //
}
