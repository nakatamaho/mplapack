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

void Cla_syrfsx_extended(INTEGER const &prec_type, const char *uplo, INTEGER const &n, INTEGER const &nrhs, COMPLEX *a, INTEGER const &lda, COMPLEX *af, INTEGER const &ldaf, INTEGER *ipiv, bool const &colequ, REAL *c, COMPLEX *b, INTEGER const &ldb, COMPLEX *y, INTEGER const &ldy, REAL *berr_out, INTEGER const &n_norms, REAL *err_bnds_norm, REAL *err_bnds_comp, COMPLEX *res, REAL *ayb, COMPLEX *dy, COMPLEX *y_tail, REAL const &rcond, INTEGER const &ithresh, REAL const &rthresh, REAL const &dz_ub, bool const &ignore_cwise, INTEGER &info) {
    COMPLEX zdum = 0.0;
    bool upper = false;
    REAL eps = 0.0;
    REAL hugeval = 0.0;
    REAL incr_thresh = 0.0;
    INTEGER uplo2 = 0;
    INTEGER j = 0;
    const INTEGER extra_residual = 1;
    INTEGER y_prec_state = 0;
    const INTEGER extra_y = 2;
    INTEGER i = 0;
    REAL dxrat = 0.0;
    REAL dxratmax = 0.0;
    REAL dzrat = 0.0;
    REAL dzratmax = 0.0;
    REAL final_dx_x = 0.0;
    REAL final_dz_z = 0.0;
    REAL prevnormdx = 0.0;
    REAL prev_dz_z = 0.0;
    REAL dz_z = 0.0;
    REAL dx_x = 0.0;
    const INTEGER working_state = 1;
    INTEGER x_state = 0;
    const INTEGER unstable_state = 0;
    INTEGER z_state = 0;
    bool incr_prec = false;
    INTEGER cnt = 0;
    const INTEGER base_residual = 0;
    REAL normx = 0.0;
    REAL normy = 0.0;
    REAL normdx = 0.0;
    REAL ymin = 0.0;
    REAL yk = 0.0;
    REAL dyk = 0.0;
    const INTEGER noprog_state = 3;
    const INTEGER conv_state = 2;
    const INTEGER la_linrx_err_i = 2;
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
    //     .. Parameters ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function Definitions ..
    abs1[zdum - 1] = abs(zdum.real()) + abs(zdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (nrhs < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldaf < max((INTEGER)1, n)) {
        info = -8;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -13;
    } else if (ldy < max((INTEGER)1, n)) {
        info = -15;
    }
    if (info != 0) {
        Mxerbla("Cla_herfsx_extended", -info);
        return;
    }
    eps = dlamch("Epsilon");
    hugeval = dlamch("Overflow");
    //     Force HUGEVAL to Inf
    hugeval = hugeval * hugeval;
    //     Using HUGEVAL may lead to spurious underflows.
    incr_thresh = n.real() * eps;
    //
    if (Mlsame(uplo, "L")) {
        uplo2 = iMlauplo["L" - 1];
    } else {
        uplo2 = iMlauplo["U" - 1];
    }
    //
    for (j = 1; j <= nrhs; j = j + 1) {
        y_prec_state = extra_residual;
        if (y_prec_state == extra_y) {
            for (i = 1; i <= n; i = i + 1) {
                y_tail[i - 1] = 0.0;
            }
        }
        //
        dxrat = 0.0;
        dxratmax = 0.0;
        dzrat = 0.0;
        dzratmax = 0.0;
        final_dx_x = hugeval;
        final_dz_z = hugeval;
        prevnormdx = hugeval;
        prev_dz_z = hugeval;
        dz_z = hugeval;
        dx_x = hugeval;
        //
        x_state = working_state;
        z_state = unstable_state;
        incr_prec = false;
        //
        for (cnt = 1; cnt <= ithresh; cnt = cnt + 1) {
            //
            //         Compute residual RES = B_s - op(A_s) * Y,
            //             op(A) = A, A**T, or A**H depending on TRANS (and type).
            //
            Ccopy(n, b[(j - 1) * ldb], 1, res, 1);
            if (y_prec_state == base_residual) {
                Csymv(uplo, n, COMPLEX(-1.0), a, lda, y[(j - 1) * ldy], 1, COMPLEX(1.0), res, 1);
            } else if (y_prec_state == extra_residual) {
                blas_Csymv_x(uplo2, n, COMPLEX(-1.0), a, lda, y[(j - 1) * ldy], 1, COMPLEX(1.0), res, 1, prec_type);
            } else {
                blas_Csymv2_x(uplo2, n, COMPLEX(-1.0), a, lda, y[(j - 1) * ldy], y_tail, 1, COMPLEX(1.0), res, 1, prec_type);
            }
            //
            //         XXX: RES is no longer needed.
            Ccopy(n, res, 1, dy, 1);
            Csytrs(uplo, n, 1, af, ldaf, ipiv, dy, n, info);
            //
            //         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.
            //
            normx = 0.0;
            normy = 0.0;
            normdx = 0.0;
            dz_z = 0.0;
            ymin = hugeval;
            //
            for (i = 1; i <= n; i = i + 1) {
                yk = abs1[y[(i - 1) + (j - 1) * ldy] - 1];
                dyk = abs1[dy[i - 1] - 1];
                //
                if (yk != 0.0) {
                    dz_z = max(dz_z, dyk / yk);
                } else if (dyk != 0.0) {
                    dz_z = hugeval;
                }
                //
                ymin = min(ymin, yk);
                //
                normy = max(normy, yk);
                //
                if (colequ) {
                    normx = max(normx, yk * c[i - 1]);
                    normdx = max(normdx, dyk * c[i - 1]);
                } else {
                    normx = normy;
                    normdx = max(normdx, dyk);
                }
            }
            //
            if (normx != 0.0) {
                dx_x = normdx / normx;
            } else if (normdx == 0.0) {
                dx_x = 0.0;
            } else {
                dx_x = hugeval;
            }
            //
            dxrat = normdx / prevnormdx;
            dzrat = dz_z / prev_dz_z;
            //
            //         Check termination criteria.
            //
            if (ymin * rcond < incr_thresh * normy && y_prec_state < extra_y) {
                incr_prec = true;
            }
            //
            if (x_state == noprog_state && dxrat <= rthresh) {
                x_state = working_state;
            }
            if (x_state == working_state) {
                if (dx_x <= eps) {
                    x_state = conv_state;
                } else if (dxrat > rthresh) {
                    if (y_prec_state != extra_y) {
                        incr_prec = true;
                    } else {
                        x_state = noprog_state;
                    }
                } else {
                    if (dxrat > dxratmax) {
                        dxratmax = dxrat;
                    }
                }
                if (x_state > working_state) {
                    final_dx_x = dx_x;
                }
            }
            //
            if (z_state == unstable_state && dz_z <= dz_ub) {
                z_state = working_state;
            }
            if (z_state == noprog_state && dzrat <= rthresh) {
                z_state = working_state;
            }
            if (z_state == working_state) {
                if (dz_z <= eps) {
                    z_state = conv_state;
                } else if (dz_z > dz_ub) {
                    z_state = unstable_state;
                    dzratmax = 0.0;
                    final_dz_z = hugeval;
                } else if (dzrat > rthresh) {
                    if (y_prec_state != extra_y) {
                        incr_prec = true;
                    } else {
                        z_state = noprog_state;
                    }
                } else {
                    if (dzrat > dzratmax) {
                        dzratmax = dzrat;
                    }
                }
                if (z_state > working_state) {
                    final_dz_z = dz_z;
                }
            }
            //
            if (x_state != working_state && (ignore_cwise || z_state != working_state)) {
                goto statement_666;
            }
            //
            if (incr_prec) {
                incr_prec = false;
                y_prec_state++;
                for (i = 1; i <= n; i = i + 1) {
                    y_tail[i - 1] = 0.0;
                }
            }
            //
            prevnormdx = normdx;
            prev_dz_z = dz_z;
            //
            //           Update soluton.
            //
            if (y_prec_state < extra_y) {
                Caxpy(n, COMPLEX(1.0), dy, 1, y[(j - 1) * ldy], 1);
            } else {
                Cla_wwaddw(n, y[(j - 1) * ldy], y_tail, dy);
            }
            //
        }
    //        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT.
    statement_666:
        //
        //     Set final_* when cnt hits ithresh.
        //
        if (x_state == working_state) {
            final_dx_x = dx_x;
        }
        if (z_state == working_state) {
            final_dz_z = dz_z;
        }
        //
        //     Compute error bounds.
        //
        if (n_norms >= 1) {
            err_bnds_norm[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_norm] = final_dx_x / (1 - dxratmax);
        }
        if (n_norms >= 2) {
            err_bnds_comp[(j - 1) + (la_linrx_err_i - 1) * lderr_bnds_comp] = final_dz_z / (1 - dzratmax);
        }
        //
        //     Compute componentwise relative backward error from formula
        //         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
        //     where abs(Z) is the componentwise absolute value of the matrix
        //     or vector Z.
        //
        //        Compute residual RES = B_s - op(A_s) * Y,
        //            op(A) = A, A**T, or A**H depending on TRANS (and type).
        //
        Ccopy(n, b[(j - 1) * ldb], 1, res, 1);
        Csymv(uplo, n, COMPLEX(-1.0), a, lda, y[(j - 1) * ldy], 1, COMPLEX(1.0), res, 1);
        //
        for (i = 1; i <= n; i = i + 1) {
            ayb[i - 1] = abs1[b[(i - 1) + (j - 1) * ldb] - 1];
        }
        //
        //     Compute abs(op(A_s))*abs(Y) + abs(B_s).
        //
        Cla_syamv(uplo2, n, 1.0, a, lda, y[(j - 1) * ldy], 1, 1.0, ayb, 1);
        //
        Cla_lin_berr(n, n, 1, res, ayb, berr_out[j - 1]);
        //
        //     End of loop for each RHS.
        //
    }
    //
}
