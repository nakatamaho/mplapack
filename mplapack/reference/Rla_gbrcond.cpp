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

REAL Rla_gbrcond(const char *trans, INTEGER const &n, INTEGER const &kl, INTEGER const &ku, REAL *ab, INTEGER const &ldab, REAL *afb, INTEGER const &ldafb, INTEGER *ipiv, INTEGER const &cmode, REAL *c, INTEGER &info, REAL *work, INTEGER *iwork) {
    REAL return_value = 0.0;
    bool notrans = false;
    INTEGER kd = 0;
    INTEGER ke = 0;
    INTEGER i = 0;
    REAL tmp = 0.0;
    INTEGER j = 0;
    REAL ainvnm = 0.0;
    INTEGER kase = 0;
    arr_1d<3, INTEGER> isave(fill0);
    //
    //  -- LAPACK computational routine --
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
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    return_value = 0.0;
    //
    info = 0;
    notrans = Mlsame(trans, "N");
    if (!notrans && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0 || kl > n - 1) {
        info = -3;
    } else if (ku < 0 || ku > n - 1) {
        info = -4;
    } else if (ldab < kl + ku + 1) {
        info = -6;
    } else if (ldafb < 2 * kl + ku + 1) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Rla_gbrcond", -info);
        return return_value;
    }
    if (n == 0) {
        return_value = 1.0;
        return return_value;
    }
    //
    //     Compute the equilibration matrix R such that
    //     inv(R)*A*C has unit 1-norm.
    //
    kd = ku + 1;
    ke = kl + 1;
    if (notrans) {
        for (i = 1; i <= n; i = i + 1) {
            tmp = 0.0;
            if (cmode == 1) {
                for (j = max(i - kl, 1); j <= min(i + ku, n); j = j + 1) {
                    tmp += abs(ab[((kd + i - j) - 1) + (j - 1) * ldab] * c[j - 1]);
                }
            } else if (cmode == 0) {
                for (j = max(i - kl, 1); j <= min(i + ku, n); j = j + 1) {
                    tmp += abs(ab[((kd + i - j) - 1) + (j - 1) * ldab]);
                }
            } else {
                for (j = max(i - kl, 1); j <= min(i + ku, n); j = j + 1) {
                    tmp += abs(ab[((kd + i - j) - 1) + (j - 1) * ldab] / c[j - 1]);
                }
            }
            work[(2 * n + i) - 1] = tmp;
        }
    } else {
        for (i = 1; i <= n; i = i + 1) {
            tmp = 0.0;
            if (cmode == 1) {
                for (j = max(i - kl, 1); j <= min(i + ku, n); j = j + 1) {
                    tmp += abs(ab[((ke - i + j) - 1) + (i - 1) * ldab] * c[j - 1]);
                }
            } else if (cmode == 0) {
                for (j = max(i - kl, 1); j <= min(i + ku, n); j = j + 1) {
                    tmp += abs(ab[((ke - i + j) - 1) + (i - 1) * ldab]);
                }
            } else {
                for (j = max(i - kl, 1); j <= min(i + ku, n); j = j + 1) {
                    tmp += abs(ab[((ke - i + j) - 1) + (i - 1) * ldab] / c[j - 1]);
                }
            }
            work[(2 * n + i) - 1] = tmp;
        }
    }
    //
    //     Estimate the norm of inv(op(A)).
    //
    ainvnm = 0.0;
    //
    kase = 0;
statement_10:
    Rlacn2(n, work[(n + 1) - 1], work, iwork, ainvnm, kase, isave);
    if (kase != 0) {
        if (kase == 2) {
            //
            //           Multiply by R.
            //
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = work[i - 1] * work[(2 * n + i) - 1];
            }
            //
            if (notrans) {
                Rgbtrs("No transpose", n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
            } else {
                Rgbtrs("Transpose", n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
            }
            //
            //           Multiply by inv(C).
            //
            if (cmode == 1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = work[i - 1] / c[i - 1];
                }
            } else if (cmode == -1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = work[i - 1] * c[i - 1];
                }
            }
        } else {
            //
            //           Multiply by inv(C**T).
            //
            if (cmode == 1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = work[i - 1] / c[i - 1];
                }
            } else if (cmode == -1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = work[i - 1] * c[i - 1];
                }
            }
            //
            if (notrans) {
                Rgbtrs("Transpose", n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
            } else {
                Rgbtrs("No transpose", n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
            }
            //
            //           Multiply by R.
            //
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = work[i - 1] * work[(2 * n + i) - 1];
            }
        }
        goto statement_10;
    }
    //
    //     Compute the estimate of the reciprocal condition number.
    //
    if (ainvnm != 0.0) {
        return_value = (1.0 / ainvnm);
    }
    //
    return return_value;
    //
}
