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

void Rlaed6(INTEGER const kniter, bool const orgati, REAL const rho, REAL *d, REAL *z, REAL const finit, REAL &tau, INTEGER &info) {
    REAL lbd = 0.0;
    REAL ubd = 0.0;
    const REAL zero = 0.0;
    INTEGER niter = 0;
    const REAL two = 2.0;
    REAL temp = 0.0;
    REAL c = 0.0;
    REAL a = 0.0;
    REAL b = 0.0;
    const REAL four = 4.0;
    REAL eps = 0.0;
    REAL base = 0.0;
    const REAL three = 3.0;
    REAL small1 = 0.0;
    const REAL one = 1.0;
    REAL sminv1 = 0.0;
    REAL small2 = 0.0;
    REAL sminv2 = 0.0;
    bool scale = false;
    REAL sclfac = 0.0;
    REAL sclinv = 0.0;
    INTEGER i = 0;
    REAL rscale[3];
    REAL cscale[3];
    REAL fc = 0.0;
    REAL df = 0.0;
    REAL ddf = 0.0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    REAL temp3 = 0.0;
    REAL f = 0.0;
    INTEGER iter = 0;
    const INTEGER maxit = 40;
    REAL eta = 0.0;
    REAL erretm = 0.0;
    REAL temp4 = 0.0;
    const REAL eight = 8.0;
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
    //     .. Parameters ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    //
    if (orgati) {
        lbd = d[2 - 1];
        ubd = d[3 - 1];
    } else {
        lbd = d[1 - 1];
        ubd = d[2 - 1];
    }
    if (finit < zero) {
        lbd = zero;
    } else {
        ubd = zero;
    }
    //
    niter = 1;
    tau = zero;
    if (kniter == 2) {
        if (orgati) {
            temp = (d[3 - 1] - d[2 - 1]) / two;
            c = rho + z[1 - 1] / ((d[1 - 1] - d[2 - 1]) - temp);
            a = c * (d[2 - 1] + d[3 - 1]) + z[2 - 1] + z[3 - 1];
            b = c * d[2 - 1] * d[3 - 1] + z[2 - 1] * d[3 - 1] + z[3 - 1] * d[2 - 1];
        } else {
            temp = (d[1 - 1] - d[2 - 1]) / two;
            c = rho + z[3 - 1] / ((d[3 - 1] - d[2 - 1]) - temp);
            a = c * (d[1 - 1] + d[2 - 1]) + z[1 - 1] + z[2 - 1];
            b = c * d[1 - 1] * d[2 - 1] + z[1 - 1] * d[2 - 1] + z[2 - 1] * d[1 - 1];
        }
        temp = max({abs(a), abs(b), abs(c)});
        a = a / temp;
        b = b / temp;
        c = c / temp;
        if (c == zero) {
            tau = b / a;
        } else if (a <= zero) {
            tau = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
        } else {
            tau = two * b / (a + sqrt(abs(a * a - four * b * c)));
        }
        if (tau < lbd || tau > ubd) {
            tau = (lbd + ubd) / two;
        }
        if (d[1 - 1] == tau || d[2 - 1] == tau || d[3 - 1] == tau) {
            tau = zero;
        } else {
            temp = finit + tau * z[1 - 1] / (d[1 - 1] * (d[1 - 1] - tau)) + tau * z[2 - 1] / (d[2 - 1] * (d[2 - 1] - tau)) + tau * z[3 - 1] / (d[3 - 1] * (d[3 - 1] - tau));
            if (temp <= zero) {
                lbd = tau;
            } else {
                ubd = tau;
            }
            if (abs(finit) <= abs(temp)) {
                tau = zero;
            }
        }
    }
    //
    //     get machine parameters for possible scaling to avoid overflow
    //
    //     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
    //     SMINV2, EPS are not SAVEd anymore between one call to the
    //     others but recomputed at each call
    //
    eps = Rlamch("Epsilon");
    base = Rlamch("Base");
    small1 = pow(base, (castINTEGER(log(Rlamch("SafMin")) / log(base) / three)));
    sminv1 = one / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;
    //
    //     Determine if scaling of inputs necessary to avoid overflow
    //     when computing 1/TEMP**3
    //
    if (orgati) {
        temp = min(abs(d[2 - 1] - tau), abs(d[3 - 1] - tau));
    } else {
        temp = min(abs(d[1 - 1] - tau), abs(d[2 - 1] - tau));
    }
    scale = false;
    if (temp <= small1) {
        scale = true;
        if (temp <= small2) {
            //
            //        Scale up by power of radix nearest 1/SAFMIN**(2/3)
            //
            sclfac = sminv2;
            sclinv = small2;
        } else {
            //
            //        Scale up by power of radix nearest 1/SAFMIN**(1/3)
            //
            sclfac = sminv1;
            sclinv = small1;
        }
        //
        //        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
        //
        for (i = 1; i <= 3; i = i + 1) {
            rscale[i - 1] = d[i - 1] * sclfac;
            cscale[i - 1] = z[i - 1] * sclfac;
        }
        tau = tau * sclfac;
        lbd = lbd * sclfac;
        ubd = ubd * sclfac;
    } else {
        //
        //        Copy D and Z to RscalE and CscalE
        //
        for (i = 1; i <= 3; i = i + 1) {
            rscale[i - 1] = d[i - 1];
            cscale[i - 1] = z[i - 1];
        }
    }
    //
    fc = zero;
    df = zero;
    ddf = zero;
    for (i = 1; i <= 3; i = i + 1) {
        temp = one / (rscale[i - 1] - tau);
        temp1 = cscale[i - 1] * temp;
        temp2 = temp1 * temp;
        temp3 = temp2 * temp;
        fc += temp1 / rscale[i - 1];
        df += temp2;
        ddf += temp3;
    }
    f = finit + tau * fc;
    //
    if (abs(f) <= zero) {
        goto statement_60;
    }
    if (f <= zero) {
        lbd = tau;
    } else {
        ubd = tau;
    }
    //
    //        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
    //                            scheme
    //
    //     It is not hard to see that
    //
    //           1) Iterations will go up monotonically
    //              if FINIT < 0;
    //
    //           2) Iterations will go down monotonically
    //              if FINIT > 0.
    //
    iter = niter + 1;
    //
    for (niter = iter; niter <= maxit; niter = niter + 1) {
        //
        if (orgati) {
            temp1 = rscale[2 - 1] - tau;
            temp2 = rscale[3 - 1] - tau;
        } else {
            temp1 = rscale[1 - 1] - tau;
            temp2 = rscale[2 - 1] - tau;
        }
        a = (temp1 + temp2) * f - temp1 * temp2 * df;
        b = temp1 * temp2 * f;
        c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
        temp = max({abs(a), abs(b), abs(c)});
        a = a / temp;
        b = b / temp;
        c = c / temp;
        if (c == zero) {
            eta = b / a;
        } else if (a <= zero) {
            eta = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
        } else {
            eta = two * b / (a + sqrt(abs(a * a - four * b * c)));
        }
        if (f * eta >= zero) {
            eta = -f / df;
        }
        //
        tau += eta;
        if (tau < lbd || tau > ubd) {
            tau = (lbd + ubd) / two;
        }
        //
        fc = zero;
        erretm = zero;
        df = zero;
        ddf = zero;
        for (i = 1; i <= 3; i = i + 1) {
            if ((rscale[i - 1] - tau) != zero) {
                temp = one / (rscale[i - 1] - tau);
                temp1 = cscale[i - 1] * temp;
                temp2 = temp1 * temp;
                temp3 = temp2 * temp;
                temp4 = temp1 / rscale[i - 1];
                fc += temp4;
                erretm += abs(temp4);
                df += temp2;
                ddf += temp3;
            } else {
                goto statement_60;
            }
        }
        f = finit + tau * fc;
        erretm = eight * (abs(finit) + abs(tau) * erretm) + abs(tau) * df;
        if ((abs(f) <= four * eps * erretm) || ((ubd - lbd) <= four * eps * abs(tau))) {
            goto statement_60;
        }
        if (f <= zero) {
            lbd = tau;
        } else {
            ubd = tau;
        }
    }
    info = 1;
statement_60:
    //
    //     Undo scaling
    //
    if (scale) {
        tau = tau * sclinv;
    }
    //
    //     End of Rlaed6
    //
}
