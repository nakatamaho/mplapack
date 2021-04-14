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

void Rlasd4(INTEGER const n, INTEGER const i, REAL *d, REAL *z, REAL *delta, REAL const rho, REAL &sigma, REAL *work, INTEGER &info) {
    const REAL one = 1.0;
    REAL eps = 0.0;
    REAL rhoinv = 0.0;
    const REAL zero = 0.0;
    REAL tau2 = 0.0;
    INTEGER ii = 0;
    INTEGER niter = 0;
    const REAL two = 2.0e+0;
    REAL temp = 0.0;
    REAL temp1 = 0.0;
    INTEGER j = 0;
    REAL psi = 0.0;
    REAL c = 0.0;
    REAL w = 0.0;
    REAL tau = 0.0;
    REAL delsq = 0.0;
    REAL a = 0.0;
    REAL b = 0.0;
    const REAL four = 4.0e+0;
    REAL dpsi = 0.0;
    REAL erretm = 0.0;
    REAL phi = 0.0;
    REAL dphi = 0.0;
    const REAL eight = 8.0e+0;
    REAL dtnsq1 = 0.0;
    REAL dtnsq = 0.0;
    REAL eta = 0.0;
    INTEGER iter = 0;
    const INTEGER maxit = 400;
    INTEGER ip1 = 0;
    REAL delsq2 = 0.0;
    REAL sq2 = 0.0;
    bool geomavg = false;
    bool orgati = false;
    REAL sglb = 0.0;
    REAL sgub = 0.0;
    const REAL ten = 10.0;
    INTEGER iim1 = 0;
    INTEGER iip1 = 0;
    bool swtch3 = false;
    REAL dw = 0.0;
    const REAL three = 3.0e+0;
    REAL dtipsq = 0.0;
    REAL dtisq = 0.0;
    REAL dtiim = 0.0;
    REAL dtiip = 0.0;
    REAL zz[3];
    REAL dd[3];
    REAL prew = 0.0;
    bool swtch = false;
    REAL temp2 = 0.0;
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Since this routine is called in an inner loop, we do no argument
    //     checking.
    //
    //     Quick return for N=1 and 2.
    //
    info = 0;
    if (n == 1) {
        //
        //        Presumably, I=1 upon entry
        //
        sigma = sqrt(d[1 - 1] * d[1 - 1] + rho * z[1 - 1] * z[1 - 1]);
        delta[1 - 1] = one;
        work[1 - 1] = one;
        return;
    }
    if (n == 2) {
        Rlasd5(i, d, z, delta, rho, sigma, work);
        return;
    }
    //
    //     Compute machine epsilon
    //
    eps = Rlamch("Epsilon");
    rhoinv = one / rho;
    tau2 = zero;
    //
    //     The case I = N
    //
    if (i == n) {
        //
        //        Initialize some basic variables
        //
        ii = n - 1;
        niter = 1;
        //
        //        Calculate initial guess
        //
        temp = rho / two;
        //
        //        If ||Z||_2 is not one, then TEMP should be set to
        //        RHO * ||Z||_2^2 / TWO
        //
        temp1 = temp / (d[n - 1] + sqrt(d[n - 1] * d[n - 1] + temp));
        for (j = 1; j <= n; j = j + 1) {
            work[j - 1] = d[j - 1] + d[n - 1] + temp1;
            delta[j - 1] = (d[j - 1] - d[n - 1]) - temp1;
        }
        //
        psi = zero;
        for (j = 1; j <= n - 2; j = j + 1) {
            psi += z[j - 1] * z[j - 1] / (delta[j - 1] * work[j - 1]);
        }
        //
        c = rhoinv + psi;
        w = c + z[ii - 1] * z[ii - 1] / (delta[ii - 1] * work[ii - 1]) + z[n - 1] * z[n - 1] / (delta[n - 1] * work[n - 1]);
        //
        if (w <= zero) {
            temp1 = sqrt(d[n - 1] * d[n - 1] + rho);
            temp = z[(n - 1) - 1] * z[(n - 1) - 1] / ((d[(n - 1) - 1] + temp1) * (d[n - 1] - d[(n - 1) - 1] + rho / (d[n - 1] + temp1))) + z[n - 1] * z[n - 1] / rho;
            //
            //           The following TAU2 is to approximate
            //           SIGMA_n^2 - D( N )*D( N )
            //
            if (c <= temp) {
                tau = rho;
            } else {
                delsq = (d[n - 1] - d[(n - 1) - 1]) * (d[n - 1] + d[(n - 1) - 1]);
                a = -c * delsq + z[(n - 1) - 1] * z[(n - 1) - 1] + z[n - 1] * z[n - 1];
                b = z[n - 1] * z[n - 1] * delsq;
                if (a < zero) {
                    tau2 = two * b / (sqrt(a * a + four * b * c) - a);
                } else {
                    tau2 = (a + sqrt(a * a + four * b * c)) / (two * c);
                }
                tau = tau2 / (d[n - 1] + sqrt(d[n - 1] * d[n - 1] + tau2));
            }
            //
            //           It can be proved that
            //               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO
            //
        } else {
            delsq = (d[n - 1] - d[(n - 1) - 1]) * (d[n - 1] + d[(n - 1) - 1]);
            a = -c * delsq + z[(n - 1) - 1] * z[(n - 1) - 1] + z[n - 1] * z[n - 1];
            b = z[n - 1] * z[n - 1] * delsq;
            //
            //           The following TAU2 is to approximate
            //           SIGMA_n^2 - D( N )*D( N )
            //
            if (a < zero) {
                tau2 = two * b / (sqrt(a * a + four * b * c) - a);
            } else {
                tau2 = (a + sqrt(a * a + four * b * c)) / (two * c);
            }
            tau = tau2 / (d[n - 1] + sqrt(d[n - 1] * d[n - 1] + tau2));
            //
            //           It can be proved that
            //           D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2
            //
        }
        //
        //        The following TAU is to approximate SIGMA_n - D( N )
        //
        //         TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )
        //
        sigma = d[n - 1] + tau;
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = (d[j - 1] - d[n - 1]) - tau;
            work[j - 1] = d[j - 1] + d[n - 1] + tau;
        }
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= ii; j = j + 1) {
            temp = z[j - 1] / (delta[j - 1] * work[j - 1]);
            psi += z[j - 1] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        //
        //        Evaluate PHI and the derivative DPHI
        //
        temp = z[n - 1] / (delta[n - 1] * work[n - 1]);
        phi = z[n - 1] * temp;
        dphi = temp * temp;
        erretm = eight * (-phi - psi) + erretm - phi + rhoinv;
        //    $          + ABS( TAU2 )*( DPSI+DPHI )
        //
        w = rhoinv + phi + psi;
        //
        //        Test for convergence
        //
        if (abs(w) <= eps * erretm) {
            goto statement_240;
        }
        //
        //        Calculate the new step
        //
        niter++;
        dtnsq1 = work[(n - 1) - 1] * delta[(n - 1) - 1];
        dtnsq = work[n - 1] * delta[n - 1];
        c = w - dtnsq1 * dpsi - dtnsq * dphi;
        a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
        b = dtnsq * dtnsq1 * w;
        if (c < zero) {
            c = abs(c);
        }
        if (c == zero) {
            eta = rho - sigma * sigma;
        } else if (a >= zero) {
            eta = (a + sqrt(abs(a * a - four * b * c))) / (two * c);
        } else {
            eta = two * b / (a - sqrt(abs(a * a - four * b * c)));
        }
        //
        //        Note, eta should be positive if w is negative, and
        //        eta should be negative otherwise. However,
        //        if for some reason caused by roundoff, eta*w > 0,
        //        we simply use one Newton step instead. This way
        //        will guarantee eta*w < 0.
        //
        if (w * eta > zero) {
            eta = -w / (dpsi + dphi);
        }
        temp = eta - dtnsq;
        if (temp > rho) {
            eta = rho + dtnsq;
        }
        //
        eta = eta / (sigma + sqrt(eta + sigma * sigma));
        tau += eta;
        sigma += eta;
        //
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = delta[j - 1] - eta;
            work[j - 1] += eta;
        }
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= ii; j = j + 1) {
            temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
            psi += z[j - 1] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        //
        //        Evaluate PHI and the derivative DPHI
        //
        tau2 = work[n - 1] * delta[n - 1];
        temp = z[n - 1] / tau2;
        phi = z[n - 1] * temp;
        dphi = temp * temp;
        erretm = eight * (-phi - psi) + erretm - phi + rhoinv;
        //    $          + ABS( TAU2 )*( DPSI+DPHI )
        //
        w = rhoinv + phi + psi;
        //
        //        Main loop to update the values of the array   DELTA
        //
        iter = niter + 1;
        //
        for (niter = iter; niter <= maxit; niter = niter + 1) {
            //
            //           Test for convergence
            //
            if (abs(w) <= eps * erretm) {
                goto statement_240;
            }
            //
            //           Calculate the new step
            //
            dtnsq1 = work[(n - 1) - 1] * delta[(n - 1) - 1];
            dtnsq = work[n - 1] * delta[n - 1];
            c = w - dtnsq1 * dpsi - dtnsq * dphi;
            a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
            b = dtnsq1 * dtnsq * w;
            if (a >= zero) {
                eta = (a + sqrt(abs(a * a - four * b * c))) / (two * c);
            } else {
                eta = two * b / (a - sqrt(abs(a * a - four * b * c)));
            }
            //
            //           Note, eta should be positive if w is negative, and
            //           eta should be negative otherwise. However,
            //           if for some reason caused by roundoff, eta*w > 0,
            //           we simply use one Newton step instead. This way
            //           will guarantee eta*w < 0.
            //
            if (w * eta > zero) {
                eta = -w / (dpsi + dphi);
            }
            temp = eta - dtnsq;
            if (temp <= zero) {
                eta = eta / two;
            }
            //
            eta = eta / (sigma + sqrt(eta + sigma * sigma));
            tau += eta;
            sigma += eta;
            //
            for (j = 1; j <= n; j = j + 1) {
                delta[j - 1] = delta[j - 1] - eta;
                work[j - 1] += eta;
            }
            //
            //           Evaluate PSI and the derivative DPSI
            //
            dpsi = zero;
            psi = zero;
            erretm = zero;
            for (j = 1; j <= ii; j = j + 1) {
                temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
                psi += z[j - 1] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            //
            //           Evaluate PHI and the derivative DPHI
            //
            tau2 = work[n - 1] * delta[n - 1];
            temp = z[n - 1] / tau2;
            phi = z[n - 1] * temp;
            dphi = temp * temp;
            erretm = eight * (-phi - psi) + erretm - phi + rhoinv;
            //    $             + ABS( TAU2 )*( DPSI+DPHI )
            //
            w = rhoinv + phi + psi;
        }
        //
        //        Return with INFO = 1, NITER = MAXIT and not converged
        //
        info = 1;
        goto statement_240;
        //
        //        End for the case I = N
        //
    } else {
        //
        //        The case for I < N
        //
        niter = 1;
        ip1 = i + 1;
        //
        //        Calculate initial guess
        //
        delsq = (d[ip1 - 1] - d[i - 1]) * (d[ip1 - 1] + d[i - 1]);
        delsq2 = delsq / two;
        sq2 = sqrt((d[i - 1] * d[i - 1] + d[ip1 - 1] * d[ip1 - 1]) / two);
        temp = delsq2 / (d[i - 1] + sq2);
        for (j = 1; j <= n; j = j + 1) {
            work[j - 1] = d[j - 1] + d[i - 1] + temp;
            delta[j - 1] = (d[j - 1] - d[i - 1]) - temp;
        }
        //
        psi = zero;
        for (j = 1; j <= i - 1; j = j + 1) {
            psi += z[j - 1] * z[j - 1] / (work[j - 1] * delta[j - 1]);
        }
        //
        phi = zero;
        for (j = n; j >= i + 2; j = j - 1) {
            phi += z[j - 1] * z[j - 1] / (work[j - 1] * delta[j - 1]);
        }
        c = rhoinv + psi + phi;
        w = c + z[i - 1] * z[i - 1] / (work[i - 1] * delta[i - 1]) + z[ip1 - 1] * z[ip1 - 1] / (work[ip1 - 1] * delta[ip1 - 1]);
        //
        geomavg = false;
        if (w > zero) {
            //
            //           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
            //
            //           We choose d(i) as origin.
            //
            orgati = true;
            ii = i;
            sglb = zero;
            sgub = delsq2 / (d[i - 1] + sq2);
            a = c * delsq + z[i - 1] * z[i - 1] + z[ip1 - 1] * z[ip1 - 1];
            b = z[i - 1] * z[i - 1] * delsq;
            if (a > zero) {
                tau2 = two * b / (a + sqrt(abs(a * a - four * b * c)));
            } else {
                tau2 = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
            }
            //
            //           TAU2 now is an estimation of SIGMA^2 - D( I )^2. The
            //           following, however, is the corresponding estimation of
            //           SIGMA - D( I ).
            //
            tau = tau2 / (d[i - 1] + sqrt(d[i - 1] * d[i - 1] + tau2));
            temp = sqrt(eps);
            if ((d[i - 1] <= temp * d[ip1 - 1]) && (abs(z[i - 1]) <= temp) && (d[i - 1] > zero)) {
                tau = min(ten * d[i - 1], sgub);
                geomavg = true;
            }
        } else {
            //
            //           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
            //
            //           We choose d(i+1) as origin.
            //
            orgati = false;
            ii = ip1;
            sglb = -delsq2 / (d[ii - 1] + sq2);
            sgub = zero;
            a = c * delsq - z[i - 1] * z[i - 1] - z[ip1 - 1] * z[ip1 - 1];
            b = z[ip1 - 1] * z[ip1 - 1] * delsq;
            if (a < zero) {
                tau2 = two * b / (a - sqrt(abs(a * a + four * b * c)));
            } else {
                tau2 = -(a + sqrt(abs(a * a + four * b * c))) / (two * c);
            }
            //
            //           TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The
            //           following, however, is the corresponding estimation of
            //           SIGMA - D( IP1 ).
            //
            tau = tau2 / (d[ip1 - 1] + sqrt(abs(d[ip1 - 1] * d[ip1 - 1] + tau2)));
        }
        //
        sigma = d[ii - 1] + tau;
        for (j = 1; j <= n; j = j + 1) {
            work[j - 1] = d[j - 1] + d[ii - 1] + tau;
            delta[j - 1] = (d[j - 1] - d[ii - 1]) - tau;
        }
        iim1 = ii - 1;
        iip1 = ii + 1;
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= iim1; j = j + 1) {
            temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
            psi += z[j - 1] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        //
        //        Evaluate PHI and the derivative DPHI
        //
        dphi = zero;
        phi = zero;
        for (j = n; j >= iip1; j = j - 1) {
            temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
            phi += z[j - 1] * temp;
            dphi += temp * temp;
            erretm += phi;
        }
        //
        w = rhoinv + phi + psi;
        //
        //        W is the value of the secular function with
        //        its ii-th element removed.
        //
        swtch3 = false;
        if (orgati) {
            if (w < zero) {
                swtch3 = true;
            }
        } else {
            if (w > zero) {
                swtch3 = true;
            }
        }
        if (ii == 1 || ii == n) {
            swtch3 = false;
        }
        //
        temp = z[ii - 1] / (work[ii - 1] * delta[ii - 1]);
        dw = dpsi + dphi + temp * temp;
        temp = z[ii - 1] * temp;
        w += temp;
        erretm = eight * (phi - psi) + erretm + two * rhoinv + three * abs(temp);
        //    $          + ABS( TAU2 )*DW
        //
        //        Test for convergence
        //
        if (abs(w) <= eps * erretm) {
            goto statement_240;
        }
        //
        if (w <= zero) {
            sglb = max(sglb, tau);
        } else {
            sgub = min(sgub, tau);
        }
        //
        //        Calculate the new step
        //
        niter++;
        if (!swtch3) {
            dtipsq = work[ip1 - 1] * delta[ip1 - 1];
            dtisq = work[i - 1] * delta[i - 1];
            if (orgati) {
                c = w - dtipsq * dw + delsq * pow2((z[i - 1] / dtisq));
            } else {
                c = w - dtisq * dw - delsq * pow2((z[ip1 - 1] / dtipsq));
            }
            a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
            b = dtipsq * dtisq * w;
            if (c == zero) {
                if (a == zero) {
                    if (orgati) {
                        a = z[i - 1] * z[i - 1] + dtipsq * dtipsq * (dpsi + dphi);
                    } else {
                        a = z[ip1 - 1] * z[ip1 - 1] + dtisq * dtisq * (dpsi + dphi);
                    }
                }
                eta = b / a;
            } else if (a <= zero) {
                eta = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
            } else {
                eta = two * b / (a + sqrt(abs(a * a - four * b * c)));
            }
        } else {
            //
            //           Interpolation using THREE most relevant poles
            //
            dtiim = work[iim1 - 1] * delta[iim1 - 1];
            dtiip = work[iip1 - 1] * delta[iip1 - 1];
            temp = rhoinv + psi + phi;
            if (orgati) {
                temp1 = z[iim1 - 1] / dtiim;
                temp1 = temp1 * temp1;
                c = (temp - dtiip * (dpsi + dphi)) - (d[iim1 - 1] - d[iip1 - 1]) * (d[iim1 - 1] + d[iip1 - 1]) * temp1;
                zz[1 - 1] = z[iim1 - 1] * z[iim1 - 1];
                if (dpsi < temp1) {
                    zz[3 - 1] = dtiip * dtiip * dphi;
                } else {
                    zz[3 - 1] = dtiip * dtiip * ((dpsi - temp1) + dphi);
                }
            } else {
                temp1 = z[iip1 - 1] / dtiip;
                temp1 = temp1 * temp1;
                c = (temp - dtiim * (dpsi + dphi)) - (d[iip1 - 1] - d[iim1 - 1]) * (d[iim1 - 1] + d[iip1 - 1]) * temp1;
                if (dphi < temp1) {
                    zz[1 - 1] = dtiim * dtiim * dpsi;
                } else {
                    zz[1 - 1] = dtiim * dtiim * (dpsi + (dphi - temp1));
                }
                zz[3 - 1] = z[iip1 - 1] * z[iip1 - 1];
            }
            zz[2 - 1] = z[ii - 1] * z[ii - 1];
            dd[1 - 1] = dtiim;
            dd[2 - 1] = delta[ii - 1] * work[ii - 1];
            dd[3 - 1] = dtiip;
            Rlaed6(niter, orgati, c, dd, zz, w, eta, info);
            //
            if (info != 0) {
                //
                //              If INFO is not 0, i.e., Rlaed6 failed, switch back
                //              to 2 pole interpolation.
                //
                swtch3 = false;
                info = 0;
                dtipsq = work[ip1 - 1] * delta[ip1 - 1];
                dtisq = work[i - 1] * delta[i - 1];
                if (orgati) {
                    c = w - dtipsq * dw + delsq * pow2((z[i - 1] / dtisq));
                } else {
                    c = w - dtisq * dw - delsq * pow2((z[ip1 - 1] / dtipsq));
                }
                a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                b = dtipsq * dtisq * w;
                if (c == zero) {
                    if (a == zero) {
                        if (orgati) {
                            a = z[i - 1] * z[i - 1] + dtipsq * dtipsq * (dpsi + dphi);
                        } else {
                            a = z[ip1 - 1] * z[ip1 - 1] + dtisq * dtisq * (dpsi + dphi);
                        }
                    }
                    eta = b / a;
                } else if (a <= zero) {
                    eta = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
                } else {
                    eta = two * b / (a + sqrt(abs(a * a - four * b * c)));
                }
            }
        }
        //
        //        Note, eta should be positive if w is negative, and
        //        eta should be negative otherwise. However,
        //        if for some reason caused by roundoff, eta*w > 0,
        //        we simply use one Newton step instead. This way
        //        will guarantee eta*w < 0.
        //
        if (w * eta >= zero) {
            eta = -w / dw;
        }
        //
        eta = eta / (sigma + sqrt(sigma * sigma + eta));
        temp = tau + eta;
        if (temp > sgub || temp < sglb) {
            if (w < zero) {
                eta = (sgub - tau) / two;
            } else {
                eta = (sglb - tau) / two;
            }
            if (geomavg) {
                if (w < zero) {
                    if (tau > zero) {
                        eta = sqrt(sgub * tau) - tau;
                    }
                } else {
                    if (sglb > zero) {
                        eta = sqrt(sglb * tau) - tau;
                    }
                }
            }
        }
        //
        prew = w;
        //
        tau += eta;
        sigma += eta;
        //
        for (j = 1; j <= n; j = j + 1) {
            work[j - 1] += eta;
            delta[j - 1] = delta[j - 1] - eta;
        }
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= iim1; j = j + 1) {
            temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
            psi += z[j - 1] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        //
        //        Evaluate PHI and the derivative DPHI
        //
        dphi = zero;
        phi = zero;
        for (j = n; j >= iip1; j = j - 1) {
            temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
            phi += z[j - 1] * temp;
            dphi += temp * temp;
            erretm += phi;
        }
        //
        tau2 = work[ii - 1] * delta[ii - 1];
        temp = z[ii - 1] / tau2;
        dw = dpsi + dphi + temp * temp;
        temp = z[ii - 1] * temp;
        w = rhoinv + phi + psi + temp;
        erretm = eight * (phi - psi) + erretm + two * rhoinv + three * abs(temp);
        //    $          + ABS( TAU2 )*DW
        //
        swtch = false;
        if (orgati) {
            if (-w > abs(prew) / ten) {
                swtch = true;
            }
        } else {
            if (w > abs(prew) / ten) {
                swtch = true;
            }
        }
        //
        //        Main loop to update the values of the array   DELTA and WORK
        //
        iter = niter + 1;
        //
        for (niter = iter; niter <= maxit; niter = niter + 1) {
            //
            //           Test for convergence
            //
            if (abs(w) <= eps * erretm) {
                //     $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN
                goto statement_240;
            }
            //
            if (w <= zero) {
                sglb = max(sglb, tau);
            } else {
                sgub = min(sgub, tau);
            }
            //
            //           Calculate the new step
            //
            if (!swtch3) {
                dtipsq = work[ip1 - 1] * delta[ip1 - 1];
                dtisq = work[i - 1] * delta[i - 1];
                if (!swtch) {
                    if (orgati) {
                        c = w - dtipsq * dw + delsq * pow2((z[i - 1] / dtisq));
                    } else {
                        c = w - dtisq * dw - delsq * pow2((z[ip1 - 1] / dtipsq));
                    }
                } else {
                    temp = z[ii - 1] / (work[ii - 1] * delta[ii - 1]);
                    if (orgati) {
                        dpsi += temp * temp;
                    } else {
                        dphi += temp * temp;
                    }
                    c = w - dtisq * dpsi - dtipsq * dphi;
                }
                a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                b = dtipsq * dtisq * w;
                if (c == zero) {
                    if (a == zero) {
                        if (!swtch) {
                            if (orgati) {
                                a = z[i - 1] * z[i - 1] + dtipsq * dtipsq * (dpsi + dphi);
                            } else {
                                a = z[ip1 - 1] * z[ip1 - 1] + dtisq * dtisq * (dpsi + dphi);
                            }
                        } else {
                            a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
                        }
                    }
                    eta = b / a;
                } else if (a <= zero) {
                    eta = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
                } else {
                    eta = two * b / (a + sqrt(abs(a * a - four * b * c)));
                }
            } else {
                //
                //              Interpolation using THREE most relevant poles
                //
                dtiim = work[iim1 - 1] * delta[iim1 - 1];
                dtiip = work[iip1 - 1] * delta[iip1 - 1];
                temp = rhoinv + psi + phi;
                if (swtch) {
                    c = temp - dtiim * dpsi - dtiip * dphi;
                    zz[1 - 1] = dtiim * dtiim * dpsi;
                    zz[3 - 1] = dtiip * dtiip * dphi;
                } else {
                    if (orgati) {
                        temp1 = z[iim1 - 1] / dtiim;
                        temp1 = temp1 * temp1;
                        temp2 = (d[iim1 - 1] - d[iip1 - 1]) * (d[iim1 - 1] + d[iip1 - 1]) * temp1;
                        c = temp - dtiip * (dpsi + dphi) - temp2;
                        zz[1 - 1] = z[iim1 - 1] * z[iim1 - 1];
                        if (dpsi < temp1) {
                            zz[3 - 1] = dtiip * dtiip * dphi;
                        } else {
                            zz[3 - 1] = dtiip * dtiip * ((dpsi - temp1) + dphi);
                        }
                    } else {
                        temp1 = z[iip1 - 1] / dtiip;
                        temp1 = temp1 * temp1;
                        temp2 = (d[iip1 - 1] - d[iim1 - 1]) * (d[iim1 - 1] + d[iip1 - 1]) * temp1;
                        c = temp - dtiim * (dpsi + dphi) - temp2;
                        if (dphi < temp1) {
                            zz[1 - 1] = dtiim * dtiim * dpsi;
                        } else {
                            zz[1 - 1] = dtiim * dtiim * (dpsi + (dphi - temp1));
                        }
                        zz[3 - 1] = z[iip1 - 1] * z[iip1 - 1];
                    }
                }
                dd[1 - 1] = dtiim;
                dd[2 - 1] = delta[ii - 1] * work[ii - 1];
                dd[3 - 1] = dtiip;
                Rlaed6(niter, orgati, c, dd, zz, w, eta, info);
                //
                if (info != 0) {
                    //
                    //                 If INFO is not 0, i.e., Rlaed6 failed, switch
                    //                 back to two pole interpolation
                    //
                    swtch3 = false;
                    info = 0;
                    dtipsq = work[ip1 - 1] * delta[ip1 - 1];
                    dtisq = work[i - 1] * delta[i - 1];
                    if (!swtch) {
                        if (orgati) {
                            c = w - dtipsq * dw + delsq * pow2((z[i - 1] / dtisq));
                        } else {
                            c = w - dtisq * dw - delsq * pow2((z[ip1 - 1] / dtipsq));
                        }
                    } else {
                        temp = z[ii - 1] / (work[ii - 1] * delta[ii - 1]);
                        if (orgati) {
                            dpsi += temp * temp;
                        } else {
                            dphi += temp * temp;
                        }
                        c = w - dtisq * dpsi - dtipsq * dphi;
                    }
                    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                    b = dtipsq * dtisq * w;
                    if (c == zero) {
                        if (a == zero) {
                            if (!swtch) {
                                if (orgati) {
                                    a = z[i - 1] * z[i - 1] + dtipsq * dtipsq * (dpsi + dphi);
                                } else {
                                    a = z[ip1 - 1] * z[ip1 - 1] + dtisq * dtisq * (dpsi + dphi);
                                }
                            } else {
                                a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
                            }
                        }
                        eta = b / a;
                    } else if (a <= zero) {
                        eta = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
                    } else {
                        eta = two * b / (a + sqrt(abs(a * a - four * b * c)));
                    }
                }
            }
            //
            //           Note, eta should be positive if w is negative, and
            //           eta should be negative otherwise. However,
            //           if for some reason caused by roundoff, eta*w > 0,
            //           we simply use one Newton step instead. This way
            //           will guarantee eta*w < 0.
            //
            if (w * eta >= zero) {
                eta = -w / dw;
            }
            //
            eta = eta / (sigma + sqrt(sigma * sigma + eta));
            temp = tau + eta;
            if (temp > sgub || temp < sglb) {
                if (w < zero) {
                    eta = (sgub - tau) / two;
                } else {
                    eta = (sglb - tau) / two;
                }
                if (geomavg) {
                    if (w < zero) {
                        if (tau > zero) {
                            eta = sqrt(sgub * tau) - tau;
                        }
                    } else {
                        if (sglb > zero) {
                            eta = sqrt(sglb * tau) - tau;
                        }
                    }
                }
            }
            //
            prew = w;
            //
            tau += eta;
            sigma += eta;
            //
            for (j = 1; j <= n; j = j + 1) {
                work[j - 1] += eta;
                delta[j - 1] = delta[j - 1] - eta;
            }
            //
            //           Evaluate PSI and the derivative DPSI
            //
            dpsi = zero;
            psi = zero;
            erretm = zero;
            for (j = 1; j <= iim1; j = j + 1) {
                temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
                psi += z[j - 1] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            //
            //           Evaluate PHI and the derivative DPHI
            //
            dphi = zero;
            phi = zero;
            for (j = n; j >= iip1; j = j - 1) {
                temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
                phi += z[j - 1] * temp;
                dphi += temp * temp;
                erretm += phi;
            }
            //
            tau2 = work[ii - 1] * delta[ii - 1];
            temp = z[ii - 1] / tau2;
            dw = dpsi + dphi + temp * temp;
            temp = z[ii - 1] * temp;
            w = rhoinv + phi + psi + temp;
            erretm = eight * (phi - psi) + erretm + two * rhoinv + three * abs(temp);
            //    $             + ABS( TAU2 )*DW
            //
            if (w * prew > zero && abs(w) > abs(prew) / ten) {
                swtch = !swtch;
            }
            //
        }
        //
        //        Return with INFO = 1, NITER = MAXIT and not converged
        //
        info = 1;
        //
    }
//
statement_240:;
    //
    //     End of Rlasd4
    //
}
