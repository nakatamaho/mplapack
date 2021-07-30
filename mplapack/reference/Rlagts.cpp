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

void Rlagts(INTEGER const job, INTEGER const n, REAL *a, REAL *b, REAL *c, REAL *d, INTEGER *in, REAL *y, REAL &tol, INTEGER &info) {
    REAL eps = 0.0;
    REAL sfmin = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    const REAL zero = 0.0;
    INTEGER k = 0;
    REAL temp = 0.0;
    REAL ak = 0.0;
    REAL absak = 0.0;
    REAL pert = 0.0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if ((abs(job) > 2) || (job == 0)) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    }
    if (info != 0) {
        Mxerbla("Rlagts", -info);
        return;
    }
    //
    if (n == 0) {
        return;
    }
    //
    eps = Rlamch("Epsilon");
    sfmin = Rlamch("Safe minimum");
    bignum = one / sfmin;
    //
    if (job < 0) {
        if (tol <= zero) {
            tol = abs(a[1 - 1]);
            if (n > 1) {
                tol = max({tol, REAL(abs(a[2 - 1])), REAL(abs(b[1 - 1]))});
            }
            for (k = 3; k <= n; k = k + 1) {
                tol = max({tol, REAL(abs(a[k - 1])), REAL(abs(b[(k - 1) - 1])), REAL(abs(d[(k - 2) - 1]))});
            }
            tol = tol * eps;
            if (tol == zero) {
                tol = eps;
            }
        }
    }
    //
    if (abs(job) == 1) {
        for (k = 2; k <= n; k = k + 1) {
            if (in[(k - 1) - 1] == 0) {
                y[k - 1] = y[k - 1] - c[(k - 1) - 1] * y[(k - 1) - 1];
            } else {
                temp = y[(k - 1) - 1];
                y[(k - 1) - 1] = y[k - 1];
                y[k - 1] = temp - c[(k - 1) - 1] * y[k - 1];
            }
        }
        if (job == 1) {
            for (k = n; k >= 1; k = k - 1) {
                if (k <= n - 2) {
                    temp = y[k - 1] - b[k - 1] * y[(k + 1) - 1] - d[k - 1] * y[(k + 2) - 1];
                } else if (k == n - 1) {
                    temp = y[k - 1] - b[k - 1] * y[(k + 1) - 1];
                } else {
                    temp = y[k - 1];
                }
                ak = a[k - 1];
                absak = abs(ak);
                if (absak < one) {
                    if (absak < sfmin) {
                        if (absak == zero || abs(temp) * sfmin > absak) {
                            info = k;
                            return;
                        } else {
                            temp = temp * bignum;
                            ak = ak * bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        info = k;
                        return;
                    }
                }
                y[k - 1] = temp / ak;
            }
        } else {
            for (k = n; k >= 1; k = k - 1) {
                if (k <= n - 2) {
                    temp = y[k - 1] - b[k - 1] * y[(k + 1) - 1] - d[k - 1] * y[(k + 2) - 1];
                } else if (k == n - 1) {
                    temp = y[k - 1] - b[k - 1] * y[(k + 1) - 1];
                } else {
                    temp = y[k - 1];
                }
                ak = a[k - 1];
                pert = sign(tol, ak);
            statement_40:
                absak = abs(ak);
                if (absak < one) {
                    if (absak < sfmin) {
                        if (absak == zero || abs(temp) * sfmin > absak) {
                            ak += pert;
                            pert = 2 * pert;
                            goto statement_40;
                        } else {
                            temp = temp * bignum;
                            ak = ak * bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        ak += pert;
                        pert = 2 * pert;
                        goto statement_40;
                    }
                }
                y[k - 1] = temp / ak;
            }
        }
    } else {
        //
        //        Come to here if  JOB = 2 or -2
        //
        if (job == 2) {
            for (k = 1; k <= n; k = k + 1) {
                if (k >= 3) {
                    temp = y[k - 1] - b[(k - 1) - 1] * y[(k - 1) - 1] - d[(k - 2) - 1] * y[(k - 2) - 1];
                } else if (k == 2) {
                    temp = y[k - 1] - b[(k - 1) - 1] * y[(k - 1) - 1];
                } else {
                    temp = y[k - 1];
                }
                ak = a[k - 1];
                absak = abs(ak);
                if (absak < one) {
                    if (absak < sfmin) {
                        if (absak == zero || abs(temp) * sfmin > absak) {
                            info = k;
                            return;
                        } else {
                            temp = temp * bignum;
                            ak = ak * bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        info = k;
                        return;
                    }
                }
                y[k - 1] = temp / ak;
            }
        } else {
            for (k = 1; k <= n; k = k + 1) {
                if (k >= 3) {
                    temp = y[k - 1] - b[(k - 1) - 1] * y[(k - 1) - 1] - d[(k - 2) - 1] * y[(k - 2) - 1];
                } else if (k == 2) {
                    temp = y[k - 1] - b[(k - 1) - 1] * y[(k - 1) - 1];
                } else {
                    temp = y[k - 1];
                }
                ak = a[k - 1];
                pert = sign(tol, ak);
            statement_70:
                absak = abs(ak);
                if (absak < one) {
                    if (absak < sfmin) {
                        if (absak == zero || abs(temp) * sfmin > absak) {
                            ak += pert;
                            pert = 2 * pert;
                            goto statement_70;
                        } else {
                            temp = temp * bignum;
                            ak = ak * bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        ak += pert;
                        pert = 2 * pert;
                        goto statement_70;
                    }
                }
                y[k - 1] = temp / ak;
            }
        }
        //
        for (k = n; k >= 2; k = k - 1) {
            if (in[(k - 1) - 1] == 0) {
                y[(k - 1) - 1] = y[(k - 1) - 1] - c[(k - 1) - 1] * y[k - 1];
            } else {
                temp = y[(k - 1) - 1];
                y[(k - 1) - 1] = y[k - 1];
                y[k - 1] = temp - c[(k - 1) - 1] * y[k - 1];
            }
        }
    }
    //
    //     End of Rlagts
    //
}
