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

void Cbbcsd(const char *jobu1, const char *jobu2, const char *jobv1t, const char *jobv2t, const char *trans, INTEGER const &m, INTEGER const &p, INTEGER const &q, REAL *theta, REAL *phi, COMPLEX *u1, INTEGER const &ldu1, COMPLEX *u2, INTEGER const &ldu2, COMPLEX *v1t, INTEGER const &ldv1t, COMPLEX *v2t, INTEGER const &ldv2t, REAL *b11d, REAL *b11e, REAL *b12d, REAL *b12e, REAL *b21d, REAL *b21e, REAL *b22d, REAL *b22e, REAL *rwork, INTEGER const &lrwork, INTEGER &info) {
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
    //  ===================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test input arguments
    //
    info = 0;
    bool lquery = lrwork == -1;
    bool wantu1 = Mlsame(jobu1, "Y");
    bool wantu2 = Mlsame(jobu2, "Y");
    bool wantv1t = Mlsame(jobv1t, "Y");
    bool wantv2t = Mlsame(jobv2t, "Y");
    bool colmajor = !Mlsame(trans, "T");
    //
    if (m < 0) {
        info = -6;
    } else if (p < 0 || p > m) {
        info = -7;
    } else if (q < 0 || q > m) {
        info = -8;
    } else if (q > p || q > m - p || q > m - q) {
        info = -8;
    } else if (wantu1 && ldu1 < p) {
        info = -12;
    } else if (wantu2 && ldu2 < m - p) {
        info = -14;
    } else if (wantv1t && ldv1t < q) {
        info = -16;
    } else if (wantv2t && ldv2t < m - q) {
        info = -18;
    }
    //
    //     Quick return if Q = 0
    //
    INTEGER lrworkmin = 0;
    if (info == 0 && q == 0) {
        lrworkmin = 1;
        rwork[1 - 1] = lrworkmin;
        return;
    }
    //
    //     Compute workspace
    //
    INTEGER iu1cs = 0;
    INTEGER iu1sn = 0;
    INTEGER iu2cs = 0;
    INTEGER iu2sn = 0;
    INTEGER iv1tcs = 0;
    INTEGER iv1tsn = 0;
    INTEGER iv2tcs = 0;
    INTEGER iv2tsn = 0;
    INTEGER lrworkopt = 0;
    if (info == 0) {
        iu1cs = 1;
        iu1sn = iu1cs + q;
        iu2cs = iu1sn + q;
        iu2sn = iu2cs + q;
        iv1tcs = iu2sn + q;
        iv1tsn = iv1tcs + q;
        iv2tcs = iv1tsn + q;
        iv2tsn = iv2tcs + q;
        lrworkopt = iv2tsn + q - 1;
        lrworkmin = lrworkopt;
        rwork[1 - 1] = lrworkopt;
        if (lrwork < lrworkmin && !lquery) {
            info = -28;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cbbcsd", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Get machine constants
    //
    REAL eps = dlamch("Epsilon");
    REAL unfl = dlamch("Safe minimum");
    const REAL ten = 10.0;
    const REAL hundred = 100.0;
    const REAL meighth = -0.125e0;
    REAL tolmul = max(ten, min(hundred, pow(eps, meighth)));
    REAL tol = tolmul * eps;
    const INTEGER maxitr = 6;
    REAL thresh = max(tol, maxitr * q * q * unfl);
    //
    //     Test for negligible sines or cosines
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    const REAL piover2 = 1.57079632679489661923132169163975144210e0;
    for (i = 1; i <= q; i = i + 1) {
        if (theta[i - 1] < thresh) {
            theta[i - 1] = zero;
        } else if (theta[i - 1] > piover2 - thresh) {
            theta[i - 1] = piover2;
        }
    }
    for (i = 1; i <= q - 1; i = i + 1) {
        if (phi[i - 1] < thresh) {
            phi[i - 1] = zero;
        } else if (phi[i - 1] > piover2 - thresh) {
            phi[i - 1] = piover2;
        }
    }
    //
    //     Initial deflation
    //
    INTEGER imax = q;
    while (imax > 1) {
        if (phi[(imax - 1) - 1] != zero) {
            break;
        }
        imax = imax - 1;
    }
    INTEGER imin = imax - 1;
    if (imin > 1) {
        while (phi[(imin - 1) - 1] != zero) {
            imin = imin - 1;
            if (imin <= 1) {
                break;
            }
        }
    }
    //
    //     Initialize iteration counter
    //
    INTEGER maxit = maxitr * q * q;
    INTEGER iter = 0;
    //
    //     Begin main iteration loop
    //
    REAL thetamax = 0.0;
    REAL thetamin = 0.0;
    REAL mu = 0.0;
    const REAL one = 1.0;
    REAL nu = 0.0;
    REAL sigma11 = 0.0;
    REAL dummy = 0.0;
    REAL sigma21 = 0.0;
    REAL temp = 0.0;
    REAL b11bulge = 0.0;
    REAL b21bulge = 0.0;
    REAL r = 0.0;
    REAL b12bulge = 0.0;
    REAL b22bulge = 0.0;
    REAL x1 = 0.0;
    REAL x2 = 0.0;
    REAL y1 = 0.0;
    REAL y2 = 0.0;
    bool restart11 = false;
    bool restart21 = false;
    bool restart12 = false;
    bool restart22 = false;
    const COMPLEX negonecomplex = (-1.0, 0.0);
    while (imax > 1) {
        //
        //        Compute the matrix entries
        //
        b11d[imin - 1] = cos[theta[imin - 1] - 1];
        b21d[imin - 1] = -sin[theta[imin - 1] - 1];
        for (i = imin; i <= imax - 1; i = i + 1) {
            b11e[i - 1] = -sin[theta[i - 1] - 1] * sin[phi[i - 1] - 1];
            b11d[(i + 1) - 1] = cos[(theta[(i + 1) - 1]) - 1] * cos[phi[i - 1] - 1];
            b12d[i - 1] = sin[theta[i - 1] - 1] * cos[phi[i - 1] - 1];
            b12e[i - 1] = cos[(theta[(i + 1) - 1]) - 1] * sin[phi[i - 1] - 1];
            b21e[i - 1] = -cos[theta[i - 1] - 1] * sin[phi[i - 1] - 1];
            b21d[(i + 1) - 1] = -sin[(theta[(i + 1) - 1]) - 1] * cos[phi[i - 1] - 1];
            b22d[i - 1] = cos[theta[i - 1] - 1] * cos[phi[i - 1] - 1];
            b22e[i - 1] = -sin[(theta[(i + 1) - 1]) - 1] * sin[phi[i - 1] - 1];
        }
        b12d[imax - 1] = sin[theta[imax - 1] - 1];
        b22d[imax - 1] = cos[theta[imax - 1] - 1];
        //
        //        Abort if not converging; otherwise, increment ITER
        //
        if (iter > maxit) {
            info = 0;
            for (i = 1; i <= q; i = i + 1) {
                if (phi[i - 1] != zero) {
                    info++;
                }
            }
            return;
        }
        //
        iter += imax - imin;
        //
        //        Compute shifts
        //
        thetamax = theta[imin - 1];
        thetamin = theta[imin - 1];
        for (i = imin + 1; i <= imax; i = i + 1) {
            if (theta[i - 1] > thetamax) {
                thetamax = theta[i - 1];
            }
            if (theta[i - 1] < thetamin) {
                thetamin = theta[i - 1];
            }
        }
        //
        if (thetamax > piover2 - thresh) {
            //
            //           Zero on diagonals of B11 and B22; induce deflation with a
            //           zero shift
            //
            mu = zero;
            nu = one;
            //
        } else if (thetamin < thresh) {
            //
            //           Zero on diagonals of B12 and B22; induce deflation with a
            //           zero shift
            //
            mu = one;
            nu = zero;
            //
        } else {
            //
            //           Compute shifts for B11 and B21 and use the lesser
            //
            Rlas2(b11d[(imax - 1) - 1], b11e[(imax - 1) - 1], b11d[imax - 1], sigma11, dummy);
            Rlas2(b21d[(imax - 1) - 1], b21e[(imax - 1) - 1], b21d[imax - 1], sigma21, dummy);
            //
            if (sigma11 <= sigma21) {
                mu = sigma11;
                nu = sqrt(one - pow2(mu));
                if (mu < thresh) {
                    mu = zero;
                    nu = one;
                }
            } else {
                nu = sigma21;
                mu = sqrt(1.0f - pow2(nu));
                if (nu < thresh) {
                    mu = one;
                    nu = zero;
                }
            }
        }
        //
        //        Rotate to produce bulges in B11 and B21
        //
        if (mu <= nu) {
            Rlartgs(b11d[imin - 1], b11e[imin - 1], mu, rwork[(iv1tcs + imin - 1) - 1], rwork[(iv1tsn + imin - 1) - 1]);
        } else {
            Rlartgs(b21d[imin - 1], b21e[imin - 1], nu, rwork[(iv1tcs + imin - 1) - 1], rwork[(iv1tsn + imin - 1) - 1]);
        }
        //
        temp = rwork[(iv1tcs + imin - 1) - 1] * b11d[imin - 1] + rwork[(iv1tsn + imin - 1) - 1] * b11e[imin - 1];
        b11e[imin - 1] = rwork[(iv1tcs + imin - 1) - 1] * b11e[imin - 1] - rwork[(iv1tsn + imin - 1) - 1] * b11d[imin - 1];
        b11d[imin - 1] = temp;
        b11bulge = rwork[(iv1tsn + imin - 1) - 1] * b11d[(imin + 1) - 1];
        b11d[(imin + 1) - 1] = rwork[(iv1tcs + imin - 1) - 1] * b11d[(imin + 1) - 1];
        temp = rwork[(iv1tcs + imin - 1) - 1] * b21d[imin - 1] + rwork[(iv1tsn + imin - 1) - 1] * b21e[imin - 1];
        b21e[imin - 1] = rwork[(iv1tcs + imin - 1) - 1] * b21e[imin - 1] - rwork[(iv1tsn + imin - 1) - 1] * b21d[imin - 1];
        b21d[imin - 1] = temp;
        b21bulge = rwork[(iv1tsn + imin - 1) - 1] * b21d[(imin + 1) - 1];
        b21d[(imin + 1) - 1] = rwork[(iv1tcs + imin - 1) - 1] * b21d[(imin + 1) - 1];
        //
        //        Compute THETA(IMIN)
        //
        theta[imin - 1] = atan2[((sqrt(pow2(b21d[imin - 1]) + pow2(b21bulge))) - 1) + ((sqrt(pow2(b11d[imin - 1]) + pow2(b11bulge))) - 1) * ldatan2];
        //
        //        Chase the bulges in B11(IMIN+1,IMIN) and B21(IMIN+1,IMIN)
        //
        if (pow2(b11d[imin - 1]) + pow2(b11bulge) > pow2(thresh)) {
            Rlartgp(b11bulge, b11d[imin - 1], rwork[(iu1sn + imin - 1) - 1], rwork[(iu1cs + imin - 1) - 1], r);
        } else if (mu <= nu) {
            Rlartgs(b11e[imin - 1], b11d[(imin + 1) - 1], mu, rwork[(iu1cs + imin - 1) - 1], rwork[(iu1sn + imin - 1) - 1]);
        } else {
            Rlartgs(b12d[imin - 1], b12e[imin - 1], nu, rwork[(iu1cs + imin - 1) - 1], rwork[(iu1sn + imin - 1) - 1]);
        }
        if (pow2(b21d[imin - 1]) + pow2(b21bulge) > pow2(thresh)) {
            Rlartgp(b21bulge, b21d[imin - 1], rwork[(iu2sn + imin - 1) - 1], rwork[(iu2cs + imin - 1) - 1], r);
        } else if (nu < mu) {
            Rlartgs(b21e[imin - 1], b21d[(imin + 1) - 1], nu, rwork[(iu2cs + imin - 1) - 1], rwork[(iu2sn + imin - 1) - 1]);
        } else {
            Rlartgs(b22d[imin - 1], b22e[imin - 1], mu, rwork[(iu2cs + imin - 1) - 1], rwork[(iu2sn + imin - 1) - 1]);
        }
        rwork[(iu2cs + imin - 1) - 1] = -rwork[(iu2cs + imin - 1) - 1];
        rwork[(iu2sn + imin - 1) - 1] = -rwork[(iu2sn + imin - 1) - 1];
        //
        temp = rwork[(iu1cs + imin - 1) - 1] * b11e[imin - 1] + rwork[(iu1sn + imin - 1) - 1] * b11d[(imin + 1) - 1];
        b11d[(imin + 1) - 1] = rwork[(iu1cs + imin - 1) - 1] * b11d[(imin + 1) - 1] - rwork[(iu1sn + imin - 1) - 1] * b11e[imin - 1];
        b11e[imin - 1] = temp;
        if (imax > imin + 1) {
            b11bulge = rwork[(iu1sn + imin - 1) - 1] * b11e[(imin + 1) - 1];
            b11e[(imin + 1) - 1] = rwork[(iu1cs + imin - 1) - 1] * b11e[(imin + 1) - 1];
        }
        temp = rwork[(iu1cs + imin - 1) - 1] * b12d[imin - 1] + rwork[(iu1sn + imin - 1) - 1] * b12e[imin - 1];
        b12e[imin - 1] = rwork[(iu1cs + imin - 1) - 1] * b12e[imin - 1] - rwork[(iu1sn + imin - 1) - 1] * b12d[imin - 1];
        b12d[imin - 1] = temp;
        b12bulge = rwork[(iu1sn + imin - 1) - 1] * b12d[(imin + 1) - 1];
        b12d[(imin + 1) - 1] = rwork[(iu1cs + imin - 1) - 1] * b12d[(imin + 1) - 1];
        temp = rwork[(iu2cs + imin - 1) - 1] * b21e[imin - 1] + rwork[(iu2sn + imin - 1) - 1] * b21d[(imin + 1) - 1];
        b21d[(imin + 1) - 1] = rwork[(iu2cs + imin - 1) - 1] * b21d[(imin + 1) - 1] - rwork[(iu2sn + imin - 1) - 1] * b21e[imin - 1];
        b21e[imin - 1] = temp;
        if (imax > imin + 1) {
            b21bulge = rwork[(iu2sn + imin - 1) - 1] * b21e[(imin + 1) - 1];
            b21e[(imin + 1) - 1] = rwork[(iu2cs + imin - 1) - 1] * b21e[(imin + 1) - 1];
        }
        temp = rwork[(iu2cs + imin - 1) - 1] * b22d[imin - 1] + rwork[(iu2sn + imin - 1) - 1] * b22e[imin - 1];
        b22e[imin - 1] = rwork[(iu2cs + imin - 1) - 1] * b22e[imin - 1] - rwork[(iu2sn + imin - 1) - 1] * b22d[imin - 1];
        b22d[imin - 1] = temp;
        b22bulge = rwork[(iu2sn + imin - 1) - 1] * b22d[(imin + 1) - 1];
        b22d[(imin + 1) - 1] = rwork[(iu2cs + imin - 1) - 1] * b22d[(imin + 1) - 1];
        //
        //        Inner loop: chase bulges from B11(IMIN,IMIN+2),
        //        B12(IMIN,IMIN+1), B21(IMIN,IMIN+2), and B22(IMIN,IMIN+1) to
        //        bottom-right
        //
        for (i = imin + 1; i <= imax - 1; i = i + 1) {
            //
            //           Compute PHI(I-1)
            //
            x1 = sin[(theta[(i - 1) - 1]) - 1] * b11e[(i - 1) - 1] + cos[(theta[(i - 1) - 1]) - 1] * b21e[(i - 1) - 1];
            x2 = sin[(theta[(i - 1) - 1]) - 1] * b11bulge + cos[(theta[(i - 1) - 1]) - 1] * b21bulge;
            y1 = sin[(theta[(i - 1) - 1]) - 1] * b12d[(i - 1) - 1] + cos[(theta[(i - 1) - 1]) - 1] * b22d[(i - 1) - 1];
            y2 = sin[(theta[(i - 1) - 1]) - 1] * b12bulge + cos[(theta[(i - 1) - 1]) - 1] * b22bulge;
            //
            phi[(i - 1) - 1] = atan2[((sqrt(pow2(x1) + pow2(x2))) - 1) + ((sqrt(pow2(y1) + pow2(y2))) - 1) * ldatan2];
            //
            //           Determine if there are bulges to chase or if a new direct
            //           summand has been reached
            //
            restart11 = pow2(b11e[(i - 1) - 1]) + pow2(b11bulge) <= pow2(thresh);
            restart21 = pow2(b21e[(i - 1) - 1]) + pow2(b21bulge) <= pow2(thresh);
            restart12 = pow2(b12d[(i - 1) - 1]) + pow2(b12bulge) <= pow2(thresh);
            restart22 = pow2(b22d[(i - 1) - 1]) + pow2(b22bulge) <= pow2(thresh);
            //
            //           If possible, chase bulges from B11(I-1,I+1), B12(I-1,I),
            //           B21(I-1,I+1), and B22(I-1,I). If necessary, restart bulge-
            //           chasing by applying the original shift again.
            //
            if (!restart11 && !restart21) {
                Rlartgp(x2, x1, rwork[(iv1tsn + i - 1) - 1], rwork[(iv1tcs + i - 1) - 1], r);
            } else if (!restart11 && restart21) {
                Rlartgp(b11bulge, b11e[(i - 1) - 1], rwork[(iv1tsn + i - 1) - 1], rwork[(iv1tcs + i - 1) - 1], r);
            } else if (restart11 && !restart21) {
                Rlartgp(b21bulge, b21e[(i - 1) - 1], rwork[(iv1tsn + i - 1) - 1], rwork[(iv1tcs + i - 1) - 1], r);
            } else if (mu <= nu) {
                Rlartgs(b11d[i - 1], b11e[i - 1], mu, rwork[(iv1tcs + i - 1) - 1], rwork[(iv1tsn + i - 1) - 1]);
            } else {
                Rlartgs(b21d[i - 1], b21e[i - 1], nu, rwork[(iv1tcs + i - 1) - 1], rwork[(iv1tsn + i - 1) - 1]);
            }
            rwork[(iv1tcs + i - 1) - 1] = -rwork[(iv1tcs + i - 1) - 1];
            rwork[(iv1tsn + i - 1) - 1] = -rwork[(iv1tsn + i - 1) - 1];
            if (!restart12 && !restart22) {
                Rlartgp(y2, y1, rwork[(iv2tsn + i - 1 - 1) - 1], rwork[(iv2tcs + i - 1 - 1) - 1], r);
            } else if (!restart12 && restart22) {
                Rlartgp(b12bulge, b12d[(i - 1) - 1], rwork[(iv2tsn + i - 1 - 1) - 1], rwork[(iv2tcs + i - 1 - 1) - 1], r);
            } else if (restart12 && !restart22) {
                Rlartgp(b22bulge, b22d[(i - 1) - 1], rwork[(iv2tsn + i - 1 - 1) - 1], rwork[(iv2tcs + i - 1 - 1) - 1], r);
            } else if (nu < mu) {
                Rlartgs(b12e[(i - 1) - 1], b12d[i - 1], nu, rwork[(iv2tcs + i - 1 - 1) - 1], rwork[(iv2tsn + i - 1 - 1) - 1]);
            } else {
                Rlartgs(b22e[(i - 1) - 1], b22d[i - 1], mu, rwork[(iv2tcs + i - 1 - 1) - 1], rwork[(iv2tsn + i - 1 - 1) - 1]);
            }
            //
            temp = rwork[(iv1tcs + i - 1) - 1] * b11d[i - 1] + rwork[(iv1tsn + i - 1) - 1] * b11e[i - 1];
            b11e[i - 1] = rwork[(iv1tcs + i - 1) - 1] * b11e[i - 1] - rwork[(iv1tsn + i - 1) - 1] * b11d[i - 1];
            b11d[i - 1] = temp;
            b11bulge = rwork[(iv1tsn + i - 1) - 1] * b11d[(i + 1) - 1];
            b11d[(i + 1) - 1] = rwork[(iv1tcs + i - 1) - 1] * b11d[(i + 1) - 1];
            temp = rwork[(iv1tcs + i - 1) - 1] * b21d[i - 1] + rwork[(iv1tsn + i - 1) - 1] * b21e[i - 1];
            b21e[i - 1] = rwork[(iv1tcs + i - 1) - 1] * b21e[i - 1] - rwork[(iv1tsn + i - 1) - 1] * b21d[i - 1];
            b21d[i - 1] = temp;
            b21bulge = rwork[(iv1tsn + i - 1) - 1] * b21d[(i + 1) - 1];
            b21d[(i + 1) - 1] = rwork[(iv1tcs + i - 1) - 1] * b21d[(i + 1) - 1];
            temp = rwork[(iv2tcs + i - 1 - 1) - 1] * b12e[(i - 1) - 1] + rwork[(iv2tsn + i - 1 - 1) - 1] * b12d[i - 1];
            b12d[i - 1] = rwork[(iv2tcs + i - 1 - 1) - 1] * b12d[i - 1] - rwork[(iv2tsn + i - 1 - 1) - 1] * b12e[(i - 1) - 1];
            b12e[(i - 1) - 1] = temp;
            b12bulge = rwork[(iv2tsn + i - 1 - 1) - 1] * b12e[i - 1];
            b12e[i - 1] = rwork[(iv2tcs + i - 1 - 1) - 1] * b12e[i - 1];
            temp = rwork[(iv2tcs + i - 1 - 1) - 1] * b22e[(i - 1) - 1] + rwork[(iv2tsn + i - 1 - 1) - 1] * b22d[i - 1];
            b22d[i - 1] = rwork[(iv2tcs + i - 1 - 1) - 1] * b22d[i - 1] - rwork[(iv2tsn + i - 1 - 1) - 1] * b22e[(i - 1) - 1];
            b22e[(i - 1) - 1] = temp;
            b22bulge = rwork[(iv2tsn + i - 1 - 1) - 1] * b22e[i - 1];
            b22e[i - 1] = rwork[(iv2tcs + i - 1 - 1) - 1] * b22e[i - 1];
            //
            //           Compute THETA(I)
            //
            x1 = cos[(phi[(i - 1) - 1]) - 1] * b11d[i - 1] + sin[(phi[(i - 1) - 1]) - 1] * b12e[(i - 1) - 1];
            x2 = cos[(phi[(i - 1) - 1]) - 1] * b11bulge + sin[(phi[(i - 1) - 1]) - 1] * b12bulge;
            y1 = cos[(phi[(i - 1) - 1]) - 1] * b21d[i - 1] + sin[(phi[(i - 1) - 1]) - 1] * b22e[(i - 1) - 1];
            y2 = cos[(phi[(i - 1) - 1]) - 1] * b21bulge + sin[(phi[(i - 1) - 1]) - 1] * b22bulge;
            //
            theta[i - 1] = atan2[((sqrt(pow2(y1) + pow2(y2))) - 1) + ((sqrt(pow2(x1) + pow2(x2))) - 1) * ldatan2];
            //
            //           Determine if there are bulges to chase or if a new direct
            //           summand has been reached
            //
            restart11 = pow2(b11d[i - 1]) + pow2(b11bulge) <= pow2(thresh);
            restart12 = pow2(b12e[(i - 1) - 1]) + pow2(b12bulge) <= pow2(thresh);
            restart21 = pow2(b21d[i - 1]) + pow2(b21bulge) <= pow2(thresh);
            restart22 = pow2(b22e[(i - 1) - 1]) + pow2(b22bulge) <= pow2(thresh);
            //
            //           If possible, chase bulges from B11(I+1,I), B12(I+1,I-1),
            //           B21(I+1,I), and B22(I+1,I-1). If necessary, restart bulge-
            //           chasing by applying the original shift again.
            //
            if (!restart11 && !restart12) {
                Rlartgp(x2, x1, rwork[(iu1sn + i - 1) - 1], rwork[(iu1cs + i - 1) - 1], r);
            } else if (!restart11 && restart12) {
                Rlartgp(b11bulge, b11d[i - 1], rwork[(iu1sn + i - 1) - 1], rwork[(iu1cs + i - 1) - 1], r);
            } else if (restart11 && !restart12) {
                Rlartgp(b12bulge, b12e[(i - 1) - 1], rwork[(iu1sn + i - 1) - 1], rwork[(iu1cs + i - 1) - 1], r);
            } else if (mu <= nu) {
                Rlartgs(b11e[i - 1], b11d[(i + 1) - 1], mu, rwork[(iu1cs + i - 1) - 1], rwork[(iu1sn + i - 1) - 1]);
            } else {
                Rlartgs(b12d[i - 1], b12e[i - 1], nu, rwork[(iu1cs + i - 1) - 1], rwork[(iu1sn + i - 1) - 1]);
            }
            if (!restart21 && !restart22) {
                Rlartgp(y2, y1, rwork[(iu2sn + i - 1) - 1], rwork[(iu2cs + i - 1) - 1], r);
            } else if (!restart21 && restart22) {
                Rlartgp(b21bulge, b21d[i - 1], rwork[(iu2sn + i - 1) - 1], rwork[(iu2cs + i - 1) - 1], r);
            } else if (restart21 && !restart22) {
                Rlartgp(b22bulge, b22e[(i - 1) - 1], rwork[(iu2sn + i - 1) - 1], rwork[(iu2cs + i - 1) - 1], r);
            } else if (nu < mu) {
                Rlartgs(b21e[i - 1], b21e[(i + 1) - 1], nu, rwork[(iu2cs + i - 1) - 1], rwork[(iu2sn + i - 1) - 1]);
            } else {
                Rlartgs(b22d[i - 1], b22e[i - 1], mu, rwork[(iu2cs + i - 1) - 1], rwork[(iu2sn + i - 1) - 1]);
            }
            rwork[(iu2cs + i - 1) - 1] = -rwork[(iu2cs + i - 1) - 1];
            rwork[(iu2sn + i - 1) - 1] = -rwork[(iu2sn + i - 1) - 1];
            //
            temp = rwork[(iu1cs + i - 1) - 1] * b11e[i - 1] + rwork[(iu1sn + i - 1) - 1] * b11d[(i + 1) - 1];
            b11d[(i + 1) - 1] = rwork[(iu1cs + i - 1) - 1] * b11d[(i + 1) - 1] - rwork[(iu1sn + i - 1) - 1] * b11e[i - 1];
            b11e[i - 1] = temp;
            if (i < imax - 1) {
                b11bulge = rwork[(iu1sn + i - 1) - 1] * b11e[(i + 1) - 1];
                b11e[(i + 1) - 1] = rwork[(iu1cs + i - 1) - 1] * b11e[(i + 1) - 1];
            }
            temp = rwork[(iu2cs + i - 1) - 1] * b21e[i - 1] + rwork[(iu2sn + i - 1) - 1] * b21d[(i + 1) - 1];
            b21d[(i + 1) - 1] = rwork[(iu2cs + i - 1) - 1] * b21d[(i + 1) - 1] - rwork[(iu2sn + i - 1) - 1] * b21e[i - 1];
            b21e[i - 1] = temp;
            if (i < imax - 1) {
                b21bulge = rwork[(iu2sn + i - 1) - 1] * b21e[(i + 1) - 1];
                b21e[(i + 1) - 1] = rwork[(iu2cs + i - 1) - 1] * b21e[(i + 1) - 1];
            }
            temp = rwork[(iu1cs + i - 1) - 1] * b12d[i - 1] + rwork[(iu1sn + i - 1) - 1] * b12e[i - 1];
            b12e[i - 1] = rwork[(iu1cs + i - 1) - 1] * b12e[i - 1] - rwork[(iu1sn + i - 1) - 1] * b12d[i - 1];
            b12d[i - 1] = temp;
            b12bulge = rwork[(iu1sn + i - 1) - 1] * b12d[(i + 1) - 1];
            b12d[(i + 1) - 1] = rwork[(iu1cs + i - 1) - 1] * b12d[(i + 1) - 1];
            temp = rwork[(iu2cs + i - 1) - 1] * b22d[i - 1] + rwork[(iu2sn + i - 1) - 1] * b22e[i - 1];
            b22e[i - 1] = rwork[(iu2cs + i - 1) - 1] * b22e[i - 1] - rwork[(iu2sn + i - 1) - 1] * b22d[i - 1];
            b22d[i - 1] = temp;
            b22bulge = rwork[(iu2sn + i - 1) - 1] * b22d[(i + 1) - 1];
            b22d[(i + 1) - 1] = rwork[(iu2cs + i - 1) - 1] * b22d[(i + 1) - 1];
            //
        }
        //
        //        Compute PHI(IMAX-1)
        //
        x1 = sin[(theta[(imax - 1) - 1]) - 1] * b11e[(imax - 1) - 1] + cos[(theta[(imax - 1) - 1]) - 1] * b21e[(imax - 1) - 1];
        y1 = sin[(theta[(imax - 1) - 1]) - 1] * b12d[(imax - 1) - 1] + cos[(theta[(imax - 1) - 1]) - 1] * b22d[(imax - 1) - 1];
        y2 = sin[(theta[(imax - 1) - 1]) - 1] * b12bulge + cos[(theta[(imax - 1) - 1]) - 1] * b22bulge;
        //
        phi[(imax - 1) - 1] = atan2[(abs(x1) - 1) + ((sqrt(pow2(y1) + pow2(y2))) - 1) * ldatan2];
        //
        //        Chase bulges from B12(IMAX-1,IMAX) and B22(IMAX-1,IMAX)
        //
        restart12 = pow2(b12d[(imax - 1) - 1]) + pow2(b12bulge) <= pow2(thresh);
        restart22 = pow2(b22d[(imax - 1) - 1]) + pow2(b22bulge) <= pow2(thresh);
        //
        if (!restart12 && !restart22) {
            Rlartgp(y2, y1, rwork[(iv2tsn + imax - 1 - 1) - 1], rwork[(iv2tcs + imax - 1 - 1) - 1], r);
        } else if (!restart12 && restart22) {
            Rlartgp(b12bulge, b12d[(imax - 1) - 1], rwork[(iv2tsn + imax - 1 - 1) - 1], rwork[(iv2tcs + imax - 1 - 1) - 1], r);
        } else if (restart12 && !restart22) {
            Rlartgp(b22bulge, b22d[(imax - 1) - 1], rwork[(iv2tsn + imax - 1 - 1) - 1], rwork[(iv2tcs + imax - 1 - 1) - 1], r);
        } else if (nu < mu) {
            Rlartgs(b12e[(imax - 1) - 1], b12d[imax - 1], nu, rwork[(iv2tcs + imax - 1 - 1) - 1], rwork[(iv2tsn + imax - 1 - 1) - 1]);
        } else {
            Rlartgs(b22e[(imax - 1) - 1], b22d[imax - 1], mu, rwork[(iv2tcs + imax - 1 - 1) - 1], rwork[(iv2tsn + imax - 1 - 1) - 1]);
        }
        //
        temp = rwork[(iv2tcs + imax - 1 - 1) - 1] * b12e[(imax - 1) - 1] + rwork[(iv2tsn + imax - 1 - 1) - 1] * b12d[imax - 1];
        b12d[imax - 1] = rwork[(iv2tcs + imax - 1 - 1) - 1] * b12d[imax - 1] - rwork[(iv2tsn + imax - 1 - 1) - 1] * b12e[(imax - 1) - 1];
        b12e[(imax - 1) - 1] = temp;
        temp = rwork[(iv2tcs + imax - 1 - 1) - 1] * b22e[(imax - 1) - 1] + rwork[(iv2tsn + imax - 1 - 1) - 1] * b22d[imax - 1];
        b22d[imax - 1] = rwork[(iv2tcs + imax - 1 - 1) - 1] * b22d[imax - 1] - rwork[(iv2tsn + imax - 1 - 1) - 1] * b22e[(imax - 1) - 1];
        b22e[(imax - 1) - 1] = temp;
        //
        //        Update singular vectors
        //
        if (wantu1) {
            if (colmajor) {
                Clasr("R", "V", "F", p, imax - imin + 1, rwork[(iu1cs + imin - 1) - 1], rwork[(iu1sn + imin - 1) - 1], u1[(imin - 1) * ldu1], ldu1);
            } else {
                Clasr("L", "V", "F", imax - imin + 1, p, rwork[(iu1cs + imin - 1) - 1], rwork[(iu1sn + imin - 1) - 1], u1[(imin - 1)], ldu1);
            }
        }
        if (wantu2) {
            if (colmajor) {
                Clasr("R", "V", "F", m - p, imax - imin + 1, rwork[(iu2cs + imin - 1) - 1], rwork[(iu2sn + imin - 1) - 1], u2[(imin - 1) * ldu2], ldu2);
            } else {
                Clasr("L", "V", "F", imax - imin + 1, m - p, rwork[(iu2cs + imin - 1) - 1], rwork[(iu2sn + imin - 1) - 1], u2[(imin - 1)], ldu2);
            }
        }
        if (wantv1t) {
            if (colmajor) {
                Clasr("L", "V", "F", imax - imin + 1, q, rwork[(iv1tcs + imin - 1) - 1], rwork[(iv1tsn + imin - 1) - 1], v1t[(imin - 1)], ldv1t);
            } else {
                Clasr("R", "V", "F", q, imax - imin + 1, rwork[(iv1tcs + imin - 1) - 1], rwork[(iv1tsn + imin - 1) - 1], v1t[(imin - 1) * ldv1t], ldv1t);
            }
        }
        if (wantv2t) {
            if (colmajor) {
                Clasr("L", "V", "F", imax - imin + 1, m - q, rwork[(iv2tcs + imin - 1) - 1], rwork[(iv2tsn + imin - 1) - 1], v2t[(imin - 1)], ldv2t);
            } else {
                Clasr("R", "V", "F", m - q, imax - imin + 1, rwork[(iv2tcs + imin - 1) - 1], rwork[(iv2tsn + imin - 1) - 1], v2t[(imin - 1) * ldv2t], ldv2t);
            }
        }
        //
        //        Fix signs on B11(IMAX-1,IMAX) and B21(IMAX-1,IMAX)
        //
        if (b11e[(imax - 1) - 1] + b21e[(imax - 1) - 1] > 0) {
            b11d[imax - 1] = -b11d[imax - 1];
            b21d[imax - 1] = -b21d[imax - 1];
            if (wantv1t) {
                if (colmajor) {
                    Cscal(q, negonecomplex, v1t[(imax - 1)], ldv1t);
                } else {
                    Cscal(q, negonecomplex, v1t[(imax - 1) * ldv1t], 1);
                }
            }
        }
        //
        //        Compute THETA(IMAX)
        //
        x1 = cos[(phi[(imax - 1) - 1]) - 1] * b11d[imax - 1] + sin[(phi[(imax - 1) - 1]) - 1] * b12e[(imax - 1) - 1];
        y1 = cos[(phi[(imax - 1) - 1]) - 1] * b21d[imax - 1] + sin[(phi[(imax - 1) - 1]) - 1] * b22e[(imax - 1) - 1];
        //
        theta[imax - 1] = atan2[(abs(y1) - 1) + (abs(x1) - 1) * ldatan2];
        //
        //        Fix signs on B11(IMAX,IMAX), B12(IMAX,IMAX-1), B21(IMAX,IMAX),
        //        and B22(IMAX,IMAX-1)
        //
        if (b11d[imax - 1] + b12e[(imax - 1) - 1] < 0) {
            b12d[imax - 1] = -b12d[imax - 1];
            if (wantu1) {
                if (colmajor) {
                    Cscal(p, negonecomplex, u1[(imax - 1) * ldu1], 1);
                } else {
                    Cscal(p, negonecomplex, u1[(imax - 1)], ldu1);
                }
            }
        }
        if (b21d[imax - 1] + b22e[(imax - 1) - 1] > 0) {
            b22d[imax - 1] = -b22d[imax - 1];
            if (wantu2) {
                if (colmajor) {
                    Cscal(m - p, negonecomplex, u2[(imax - 1) * ldu2], 1);
                } else {
                    Cscal(m - p, negonecomplex, u2[(imax - 1)], ldu2);
                }
            }
        }
        //
        //        Fix signs on B12(IMAX,IMAX) and B22(IMAX,IMAX)
        //
        if (b12d[imax - 1] + b22d[imax - 1] < 0) {
            if (wantv2t) {
                if (colmajor) {
                    Cscal(m - q, negonecomplex, v2t[(imax - 1)], ldv2t);
                } else {
                    Cscal(m - q, negonecomplex, v2t[(imax - 1) * ldv2t], 1);
                }
            }
        }
        //
        //        Test for negligible sines or cosines
        //
        for (i = imin; i <= imax; i = i + 1) {
            if (theta[i - 1] < thresh) {
                theta[i - 1] = zero;
            } else if (theta[i - 1] > piover2 - thresh) {
                theta[i - 1] = piover2;
            }
        }
        for (i = imin; i <= imax - 1; i = i + 1) {
            if (phi[i - 1] < thresh) {
                phi[i - 1] = zero;
            } else if (phi[i - 1] > piover2 - thresh) {
                phi[i - 1] = piover2;
            }
        }
        //
        //        Deflate
        //
        if (imax > 1) {
            while (phi[(imax - 1) - 1] == zero) {
                imax = imax - 1;
                if (imax <= 1) {
                    break;
                }
            }
        }
        if (imin > imax - 1) {
            imin = imax - 1;
        }
        if (imin > 1) {
            while (phi[(imin - 1) - 1] != zero) {
                imin = imin - 1;
                if (imin <= 1) {
                    break;
                }
            }
        }
        //
        //        Repeat main iteration loop
        //
    }
    //
    //     Postprocessing: order THETA from least to greatest
    //
    INTEGER mini = 0;
    INTEGER j = 0;
    for (i = 1; i <= q; i = i + 1) {
        //
        mini = i;
        thetamin = theta[i - 1];
        for (j = i + 1; j <= q; j = j + 1) {
            if (theta[j - 1] < thetamin) {
                mini = j;
                thetamin = theta[j - 1];
            }
        }
        //
        if (mini != i) {
            theta[mini - 1] = theta[i - 1];
            theta[i - 1] = thetamin;
            if (colmajor) {
                if (wantu1) {
                    Cswap(p, u1[(i - 1) * ldu1], 1, u1[(mini - 1) * ldu1], 1);
                }
                if (wantu2) {
                    Cswap(m - p, u2[(i - 1) * ldu2], 1, u2[(mini - 1) * ldu2], 1);
                }
                if (wantv1t) {
                    Cswap(q, v1t[(i - 1)], ldv1t, v1t[(mini - 1)], ldv1t);
                }
                if (wantv2t) {
                    Cswap(m - q, v2t[(i - 1)], ldv2t, v2t[(mini - 1)], ldv2t);
                }
            } else {
                if (wantu1) {
                    Cswap(p, u1[(i - 1)], ldu1, u1[(mini - 1)], ldu1);
                }
                if (wantu2) {
                    Cswap(m - p, u2[(i - 1)], ldu2, u2[(mini - 1)], ldu2);
                }
                if (wantv1t) {
                    Cswap(q, v1t[(i - 1) * ldv1t], 1, v1t[(mini - 1) * ldv1t], 1);
                }
                if (wantv2t) {
                    Cswap(m - q, v2t[(i - 1) * ldv2t], 1, v2t[(mini - 1) * ldv2t], 1);
                }
            }
        }
        //
    }
    //
    //     End of Cbbcsd
    //
}
