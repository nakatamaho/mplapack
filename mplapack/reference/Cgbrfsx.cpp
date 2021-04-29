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

void Cgbrfsx(const char *trans, const char *equed, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, COMPLEX *ab, INTEGER const ldab, COMPLEX *afb, INTEGER const ldafb, INTEGER *ipiv, REAL *r, REAL *c, COMPLEX *b, INTEGER const ldb, COMPLEX *x, INTEGER const ldx, REAL &rcond, REAL *berr, INTEGER const n_err_bnds, REAL *err_bnds_norm, REAL *err_bnds_comp, INTEGER const nparams, REAL *params, COMPLEX *work, REAL *rwork, INTEGER &info) {
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
    //  ==================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Check the input parameters.
    //
    INTEGER lderr_bnds_norm = nrhs;
    INTEGER lderr_bnds_comp = nrhs;  
    info = 0;
    INTEGER trans_type = iMlatrans(trans);
    const REAL itref_default = 1.0;
    INTEGER ref_type = castINTEGER(itref_default);
    const INTEGER la_linrx_itref_i = 1;
    if (nparams >= la_linrx_itref_i) {
        if (params[la_linrx_itref_i - 1] < 0.0) {
            params[la_linrx_itref_i - 1] = itref_default;
        } else {
            ref_type = params[la_linrx_itref_i - 1];
        }
    }
    //
    //     Set default parameters.
    //
    REAL illrcond_thresh = castREAL(n) * Rlamch("Epsilon");
    const REAL ithresh_default = 10.0;
    INTEGER ithresh = castINTEGER(ithresh_default);
    const REAL rthresh_default = 0.5e+0;
    REAL rthresh = rthresh_default;
    const REAL dzthresh_default = 0.25e+0;
    REAL unstable_thresh = dzthresh_default;
    const REAL componentwise_default = 1.0;
    bool ignore_cwise = componentwise_default == 0.0;
    //
    const INTEGER la_linrx_ithresh_i = 2;
    if (nparams >= la_linrx_ithresh_i) {
        if (params[la_linrx_ithresh_i - 1] < 0.0) {
            params[la_linrx_ithresh_i - 1] = ithresh;
        } else {
            ithresh = castINTEGER(params[la_linrx_ithresh_i - 1]);
        }
    }
    const INTEGER la_linrx_cwise_i = 3;
    if (nparams >= la_linrx_cwise_i) {
        if (params[la_linrx_cwise_i - 1] < 0.0) {
            if (ignore_cwise) {
                params[la_linrx_cwise_i - 1] = 0.0;
            } else {
                params[la_linrx_cwise_i - 1] = 1.0;
            }
        } else {
            ignore_cwise = params[la_linrx_cwise_i - 1] == 0.0;
        }
    }
    INTEGER n_norms = 0;
    if (ref_type == 0 || n_err_bnds == 0) {
        n_norms = 0;
    } else if (ignore_cwise) {
        n_norms = 1;
    } else {
        n_norms = 2;
    }
    //
    bool notran = Mlsame(trans, "N");
    bool rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
    bool colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
    //
    //     Test input parameters.
    //
    if (trans_type == -1) {
        info = -1;
    } else if (!rowequ && !colequ && !Mlsame(equed, "N")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (kl < 0) {
        info = -4;
    } else if (ku < 0) {
        info = -5;
    } else if (nrhs < 0) {
        info = -6;
    } else if (ldab < kl + ku + 1) {
        info = -8;
    } else if (ldafb < 2 * kl + ku + 1) {
        info = -10;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -13;
    } else if (ldx < max((INTEGER)1, n)) {
        info = -15;
    }
    if (info != 0) {
        Mxerbla("Cgbrfsx", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    INTEGER j = 0;
    const INTEGER la_linrx_trust_i = 1;
    const INTEGER la_linrx_err_i = 2;
    const INTEGER la_linrx_rcond_i = 3;
    if (n == 0 || nrhs == 0) {
        rcond = 1.0;
        for (j = 1; j <= nrhs; j = j + 1) {
            berr[j - 1] = 0.0;
            if (n_err_bnds >= 1) {
                err_bnds_norm[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_norm] = 1.0;
                err_bnds_comp[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_comp] = 1.0;
            }
            if (n_err_bnds >= 2) {
                err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] = 0.0;
                err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] = 0.0;
            }
            if (n_err_bnds >= 3) {
                err_bnds_norm[(j - 1) + (la_linrx_rcond_i - 1) * lderr_bnds_norm] = 1.0;
                err_bnds_comp[(j - 1) + (la_linrx_rcond_i - 1) * lderr_bnds_comp] = 1.0;
            }
        }
        return;
    }
    //
    //     Default to failure.
    //
    rcond = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        berr[j - 1] = 1.0;
        if (n_err_bnds >= 1) {
            err_bnds_norm[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_norm] = 1.0;
            err_bnds_comp[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_comp] = 1.0;
        }
        if (n_err_bnds >= 2) {
            err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] = 1.0;
            err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] = 1.0;
        }
        if (n_err_bnds >= 3) {
            err_bnds_norm[(j - 1) + (la_linrx_rcond_i - 1) * lderr_bnds_norm] = 0.0;
            err_bnds_comp[(j - 1) + (la_linrx_rcond_i - 1) * lderr_bnds_comp] = 0.0;
        }
    }
    //
    //     Compute the norm of A and the reciprocal of the condition
    //     number of A.
    //
    char norm ;
    if (notran) {
        norm = 'I';
    } else {
        norm = '1';
    }
    REAL anorm = Clangb(&norm, n, kl, ku, ab, ldab, rwork);
    Cgbcon(&norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, rwork, info);
    //
    //     Perform refinement on each right-hand side
    //
    INTEGER prec_type = 0;
    const REAL zero = 0.0;
    if (ref_type != 0 && info == 0) {
        //
       g prec_type = iMlaprec("E");
        //
        if (notran) {
            Cla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, rwork, &work[(n + 1) - 1], transfer(rwork[((2 * n) - 1) * ldrwork], (/ (zero, zero) /), n), rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
        } else {
            Cla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, rowequ, r, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, rwork, &work[(n + 1) - 1], transfer(rwork[((2 * n) - 1) * ldrwork], (/ (zero, zero) /), n), rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
        }
    }
    //
    REAL err_lbnd = max((INTEGER)10.0, sqrt(n.real())) * Rlamch("Epsilon");
    REAL rcond_tmp = 0.0;
    if (n_err_bnds >= 1 && n_norms >= 1) {
        //
        //     Compute scaled normwise condition number cond(A*C).
        //
        if (colequ && notran) {
            rcond_tmp = Cla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, true, info, work, rwork);
        } else if (rowequ && !notran) {
            rcond_tmp = Cla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, r, true, info, work, rwork);
        } else {
            rcond_tmp = Cla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, false, info, work, rwork);
        }
        for (j = 1; j <= nrhs; j = j + 1) {
            //
            //     Cap the error at 1.0.
            //
            if (n_err_bnds >= la_linrx_err_i && err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] > 1.0) {
                err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] = 1.0;
            }
            //
            //     Threshold the error (see LAWN).
            //
            if (rcond_tmp < illrcond_thresh) {
                err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] = 1.0;
                err_bnds_norm[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_norm] = 0.0;
                if (info <= n) {
                    info = n + j;
                }
            } else if (err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] < err_lbnd) {
                err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] = err_lbnd;
                err_bnds_norm[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_norm] = 1.0;
            }
            //
            //     Save the condition number.
            //
            if (n_err_bnds >= la_linrx_rcond_i) {
                err_bnds_norm[(j - 1) + (la_linrx_rcond_i - 1) * lderr_bnds_norm] = rcond_tmp;
            }
            //
        }
    }
    //
    REAL cwise_wrong = 0.0;
    if (n_err_bnds >= 1 && n_norms >= 2) {
        //
        //     Compute componentwise condition number cond(A*diag(Y(:,J))) for
        //     each right-hand side using the current solution as an estimate of
        //     the true solution.  If the componentwise error estimate is too
        //     large, then the solution is a lousy estimate of truth and the
        //     estimated RCOND may be too optimistic.  To avoid misleading users,
        //     the inverse condition number is set to 0.0 when the estimated
        //     cwise error is at least CWISE_WRONG.
        //
        cwise_wrong = sqrt(Rlamch("Epsilon"));
        for (j = 1; j <= nrhs; j = j + 1) {
            if (err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] < cwise_wrong) {
                rcond_tmp = Cla_gbrcond_x(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, &x[(j - 1) * ldx], info, work, rwork);
            } else {
                rcond_tmp = 0.0;
            }
            //
            //     Cap the error at 1.0.
            //
            if (n_err_bnds >= la_linrx_err_i && err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] > 1.0) {
                err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] = 1.0;
            }
            //
            //     Threshold the error (see LAWN).
            //
            if (rcond_tmp < illrcond_thresh) {
                err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] = 1.0;
                err_bnds_comp[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_comp] = 0.0;
                if (params[la_linrx_cwise_i - 1] == 1.0 && info < n + j) {
                    info = n + j;
                }
            } else if (err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] < err_lbnd) {
                err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] = err_lbnd;
                err_bnds_comp[(j - 1) + (la_linrx_trust_i - 1) * lderr_bnds_comp] = 1.0;
            }
            //
            //     Save the condition number.
            //
            if (n_err_bnds >= la_linrx_rcond_i) {
                err_bnds_comp[(j - 1) + (la_linrx_rcond_i - 1) * lderr_bnds_comp] = rcond_tmp;
            }
            //
        }
    }
    //
    //     End of Cgbrfsx
    //
}
