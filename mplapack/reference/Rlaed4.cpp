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

void Rlaed4(INTEGER const n, INTEGER const i, REAL *d, REAL *z, REAL *delta, REAL const rho, REAL &dlam, INTEGER &info) {
    const REAL one = 1.0;
    REAL eps = 0.0;
    REAL rhoinv = 0.0;
    INTEGER ii = 0;
    INTEGER niter = 0;
    const REAL two = 2.0;
    REAL midpt = 0.0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL psi = 0.0;
    REAL c = 0.0;
    REAL w = 0.0;
    REAL temp = 0.0;
    REAL tau = 0.0;
    REAL del = 0.0;
    REAL a = 0.0;
    REAL b = 0.0;
    const REAL four = 4.0;
    REAL dltlb = 0.0;
    REAL dltub = 0.0;
    REAL dpsi = 0.0;
    REAL erretm = 0.0;
    REAL phi = 0.0;
    REAL dphi = 0.0;
    const REAL eight = 8.0;
    REAL eta = 0.0;
    INTEGER iter = 0;
    const INTEGER maxit = 30;
    INTEGER ip1 = 0;
    bool orgati = false;
    INTEGER iim1 = 0;
    INTEGER iip1 = 0;
    bool swtch3 = false;
    REAL dw = 0.0;
    const REAL three = 3.0;
    REAL temp1 = 0.0;
    arr_1d<3, REAL> zz(fill0);
    REAL prew = 0.0;
    bool swtch = false;
    const REAL ten = 10.0;
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
    //     Since this routine is called in an inner loop, we do no argument
    //     checking.
    //
    //     Quick return for N=1 and 2.
    //
    info = 0;
    if (n == 1) {
        //
        //         Presumably, I=1 upon entry
        //
        dlam = d[1 - 1] + rho * z[1 - 1] * z[1 - 1];
        delta[1 - 1] = one;
        return;
    }
    if (n == 2) {
        Rlaed5(i, d, z, delta, rho, dlam);
        return;
    }
    //
    //     Compute machine epsilon
    //
    eps = Rlamch("Epsilon");
    rhoinv = one / rho;
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
        midpt = rho / two;
        //
        //        If ||Z||_2 is not one, then TEMP should be set to
        //        RHO * ||Z||_2^2 / TWO
        //
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = (d[j - 1] - d[i - 1]) - midpt;
        }
        //
        psi = zero;
        for (j = 1; j <= n - 2; j = j + 1) {
            psi += z[j - 1] * z[j - 1] / delta[j - 1];
        }
        //
        c = rhoinv + psi;
        w = c + z[ii - 1] * z[ii - 1] / delta[ii - 1] + z[n - 1] * z[n - 1] / delta[n - 1];
        //
        if (w <= zero) {
            temp = z[(n - 1) - 1] * z[(n - 1) - 1] / (d[n - 1] - d[(n - 1) - 1] + rho) + z[n - 1] * z[n - 1] / rho;
            if (c <= temp) {
                tau = rho;
            } else {
                del = d[n - 1] - d[(n - 1) - 1];
                a = -c * del + z[(n - 1) - 1] * z[(n - 1) - 1] + z[n - 1] * z[n - 1];
                b = z[n - 1] * z[n - 1] * del;
                if (a < zero) {
                    tau = two * b / (sqrt(a * a + four * b * c) - a);
                } else {
                    tau = (a + sqrt(a * a + four * b * c)) / (two * c);
                }
            }
            //
            //           It can be proved that
            //               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
            //
            dltlb = midpt;
            dltub = rho;
        } else {
            del = d[n - 1] - d[(n - 1) - 1];
            a = -c * del + z[(n - 1) - 1] * z[(n - 1) - 1] + z[n - 1] * z[n - 1];
            b = z[n - 1] * z[n - 1] * del;
            if (a < zero) {
                tau = two * b / (sqrt(a * a + four * b * c) - a);
            } else {
                tau = (a + sqrt(a * a + four * b * c)) / (two * c);
            }
            //
            //           It can be proved that
            //               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
            //
            dltlb = zero;
            dltub = midpt;
        }
        //
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = (d[j - 1] - d[i - 1]) - tau;
        }
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= ii; j = j + 1) {
            temp = z[j - 1] / delta[j - 1];
            psi += z[j - 1] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        //
        //        Evaluate PHI and the derivative DPHI
        //
        temp = z[n - 1] / delta[n - 1];
        phi = z[n - 1] * temp;
        dphi = temp * temp;
        erretm = eight * (-phi - psi) + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
        //
        w = rhoinv + phi + psi;
        //
        //        Test for convergence
        //
        if (abs(w) <= eps * erretm) {
            dlam = d[i - 1] + tau;
            goto statement_250;
        }
        //
        if (w <= zero) {
            dltlb = max(dltlb, tau);
        } else {
            dltub = min(dltub, tau);
        }
        //
        //        Calculate the new step
        //
        niter++;
        c = w - delta[(n - 1) - 1] * dpsi - delta[n - 1] * dphi;
        a = (delta[(n - 1) - 1] + delta[n - 1]) * w - delta[(n - 1) - 1] * delta[n - 1] * (dpsi + dphi);
        b = delta[(n - 1) - 1] * delta[n - 1] * w;
        if (c < zero) {
            c = abs(c);
        }
        if (c == zero) {
            //          ETA = B/A
            //           ETA = RHO - TAU
            eta = dltub - tau;
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
        temp = tau + eta;
        if (temp > dltub || temp < dltlb) {
            if (w < zero) {
                eta = (dltub - tau) / two;
            } else {
                eta = (dltlb - tau) / two;
            }
        }
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = delta[j - 1] - eta;
        }
        //
        tau += eta;
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= ii; j = j + 1) {
            temp = z[j - 1] / delta[j - 1];
            psi += z[j - 1] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        //
        //        Evaluate PHI and the derivative DPHI
        //
        temp = z[n - 1] / delta[n - 1];
        phi = z[n - 1] * temp;
        dphi = temp * temp;
        erretm = eight * (-phi - psi) + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
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
                dlam = d[i - 1] + tau;
                goto statement_250;
            }
            //
            if (w <= zero) {
                dltlb = max(dltlb, tau);
            } else {
                dltub = min(dltub, tau);
            }
            //
            //           Calculate the new step
            //
            c = w - delta[(n - 1) - 1] * dpsi - delta[n - 1] * dphi;
            a = (delta[(n - 1) - 1] + delta[n - 1]) * w - delta[(n - 1) - 1] * delta[n - 1] * (dpsi + dphi);
            b = delta[(n - 1) - 1] * delta[n - 1] * w;
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
            temp = tau + eta;
            if (temp > dltub || temp < dltlb) {
                if (w < zero) {
                    eta = (dltub - tau) / two;
                } else {
                    eta = (dltlb - tau) / two;
                }
            }
            for (j = 1; j <= n; j = j + 1) {
                delta[j - 1] = delta[j - 1] - eta;
            }
            //
            tau += eta;
            //
            //           Evaluate PSI and the derivative DPSI
            //
            dpsi = zero;
            psi = zero;
            erretm = zero;
            for (j = 1; j <= ii; j = j + 1) {
                temp = z[j - 1] / delta[j - 1];
                psi += z[j - 1] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            //
            //           Evaluate PHI and the derivative DPHI
            //
            temp = z[n - 1] / delta[n - 1];
            phi = z[n - 1] * temp;
            dphi = temp * temp;
            erretm = eight * (-phi - psi) + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
            //
            w = rhoinv + phi + psi;
        }
        //
        //        Return with INFO = 1, NITER = MAXIT and not converged
        //
        info = 1;
        dlam = d[i - 1] + tau;
        goto statement_250;
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
        del = d[ip1 - 1] - d[i - 1];
        midpt = del / two;
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = (d[j - 1] - d[i - 1]) - midpt;
        }
        //
        psi = zero;
        for (j = 1; j <= i - 1; j = j + 1) {
            psi += z[j - 1] * z[j - 1] / delta[j - 1];
        }
        //
        phi = zero;
        for (j = n; j >= i + 2; j = j - 1) {
            phi += z[j - 1] * z[j - 1] / delta[j - 1];
        }
        c = rhoinv + psi + phi;
        w = c + z[i - 1] * z[i - 1] / delta[i - 1] + z[ip1 - 1] * z[ip1 - 1] / delta[ip1 - 1];
        //
        if (w > zero) {
            //
            //           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
            //
            //           We choose d(i) as origin.
            //
            orgati = true;
            a = c * del + z[i - 1] * z[i - 1] + z[ip1 - 1] * z[ip1 - 1];
            b = z[i - 1] * z[i - 1] * del;
            if (a > zero) {
                tau = two * b / (a + sqrt(abs(a * a - four * b * c)));
            } else {
                tau = (a - sqrt(abs(a * a - four * b * c))) / (two * c);
            }
            dltlb = zero;
            dltub = midpt;
        } else {
            //
            //           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
            //
            //           We choose d(i+1) as origin.
            //
            orgati = false;
            a = c * del - z[i - 1] * z[i - 1] - z[ip1 - 1] * z[ip1 - 1];
            b = z[ip1 - 1] * z[ip1 - 1] * del;
            if (a < zero) {
                tau = two * b / (a - sqrt(abs(a * a + four * b * c)));
            } else {
                tau = -(a + sqrt(abs(a * a + four * b * c))) / (two * c);
            }
            dltlb = -midpt;
            dltub = zero;
        }
        //
        if (orgati) {
            for (j = 1; j <= n; j = j + 1) {
                delta[j - 1] = (d[j - 1] - d[i - 1]) - tau;
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                delta[j - 1] = (d[j - 1] - d[ip1 - 1]) - tau;
            }
        }
        if (orgati) {
            ii = i;
        } else {
            ii = i + 1;
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
            temp = z[j - 1] / delta[j - 1];
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
            temp = z[j - 1] / delta[j - 1];
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
        temp = z[ii - 1] / delta[ii - 1];
        dw = dpsi + dphi + temp * temp;
        temp = z[ii - 1] * temp;
        w += temp;
        erretm = eight * (phi - psi) + erretm + two * rhoinv + three * abs(temp) + abs(tau) * dw;
        //
        //        Test for convergence
        //
        if (abs(w) <= eps * erretm) {
            if (orgati) {
                dlam = d[i - 1] + tau;
            } else {
                dlam = d[ip1 - 1] + tau;
            }
            goto statement_250;
        }
        //
        if (w <= zero) {
            dltlb = max(dltlb, tau);
        } else {
            dltub = min(dltub, tau);
        }
        //
        //        Calculate the new step
        //
        niter++;
        if (!swtch3) {
            if (orgati) {
                c = w - delta[ip1 - 1] * dw - (d[i - 1] - d[ip1 - 1]) * pow2((z[i - 1] / delta[i - 1]));
            } else {
                c = w - delta[i - 1] * dw - (d[ip1 - 1] - d[i - 1]) * pow2((z[ip1 - 1] / delta[ip1 - 1]));
            }
            a = (delta[i - 1] + delta[ip1 - 1]) * w - delta[i - 1] * delta[ip1 - 1] * dw;
            b = delta[i - 1] * delta[ip1 - 1] * w;
            if (c == zero) {
                if (a == zero) {
                    if (orgati) {
                        a = z[i - 1] * z[i - 1] + delta[ip1 - 1] * delta[ip1 - 1] * (dpsi + dphi);
                    } else {
                        a = z[ip1 - 1] * z[ip1 - 1] + delta[i - 1] * delta[i - 1] * (dpsi + dphi);
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
            temp = rhoinv + psi + phi;
            if (orgati) {
                temp1 = z[iim1 - 1] / delta[iim1 - 1];
                temp1 = temp1 * temp1;
                c = temp - delta[iip1 - 1] * (dpsi + dphi) - (d[iim1 - 1] - d[iip1 - 1]) * temp1;
                zz[1 - 1] = z[iim1 - 1] * z[iim1 - 1];
                zz[3 - 1] = delta[iip1 - 1] * delta[iip1 - 1] * ((dpsi - temp1) + dphi);
            } else {
                temp1 = z[iip1 - 1] / delta[iip1 - 1];
                temp1 = temp1 * temp1;
                c = temp - delta[iim1 - 1] * (dpsi + dphi) - (d[iip1 - 1] - d[iim1 - 1]) * temp1;
                zz[1 - 1] = delta[iim1 - 1] * delta[iim1 - 1] * (dpsi + (dphi - temp1));
                zz[3 - 1] = z[iip1 - 1] * z[iip1 - 1];
            }
            zz[2 - 1] = z[ii - 1] * z[ii - 1];
            Rlaed6(niter, orgati, c, delta[iim1 - 1], zz, w, eta, info);
            if (info != 0) {
                goto statement_250;
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
        temp = tau + eta;
        if (temp > dltub || temp < dltlb) {
            if (w < zero) {
                eta = (dltub - tau) / two;
            } else {
                eta = (dltlb - tau) / two;
            }
        }
        //
        prew = w;
        //
        for (j = 1; j <= n; j = j + 1) {
            delta[j - 1] = delta[j - 1] - eta;
        }
        //
        //        Evaluate PSI and the derivative DPSI
        //
        dpsi = zero;
        psi = zero;
        erretm = zero;
        for (j = 1; j <= iim1; j = j + 1) {
            temp = z[j - 1] / delta[j - 1];
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
            temp = z[j - 1] / delta[j - 1];
            phi += z[j - 1] * temp;
            dphi += temp * temp;
            erretm += phi;
        }
        //
        temp = z[ii - 1] / delta[ii - 1];
        dw = dpsi + dphi + temp * temp;
        temp = z[ii - 1] * temp;
        w = rhoinv + phi + psi + temp;
        erretm = eight * (phi - psi) + erretm + two * rhoinv + three * abs(temp) + abs(tau + eta) * dw;
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
        tau += eta;
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
                if (orgati) {
                    dlam = d[i - 1] + tau;
                } else {
                    dlam = d[ip1 - 1] + tau;
                }
                goto statement_250;
            }
            //
            if (w <= zero) {
                dltlb = max(dltlb, tau);
            } else {
                dltub = min(dltub, tau);
            }
            //
            //           Calculate the new step
            //
            if (!swtch3) {
                if (!swtch) {
                    if (orgati) {
                        c = w - delta[ip1 - 1] * dw - (d[i - 1] - d[ip1 - 1]) * pow2((z[i - 1] / delta[i - 1]));
                    } else {
                        c = w - delta[i - 1] * dw - (d[ip1 - 1] - d[i - 1]) * pow2((z[ip1 - 1] / delta[ip1 - 1]));
                    }
                } else {
                    temp = z[ii - 1] / delta[ii - 1];
                    if (orgati) {
                        dpsi += temp * temp;
                    } else {
                        dphi += temp * temp;
                    }
                    c = w - delta[i - 1] * dpsi - delta[ip1 - 1] * dphi;
                }
                a = (delta[i - 1] + delta[ip1 - 1]) * w - delta[i - 1] * delta[ip1 - 1] * dw;
                b = delta[i - 1] * delta[ip1 - 1] * w;
                if (c == zero) {
                    if (a == zero) {
                        if (!swtch) {
                            if (orgati) {
                                a = z[i - 1] * z[i - 1] + delta[ip1 - 1] * delta[ip1 - 1] * (dpsi + dphi);
                            } else {
                                a = z[ip1 - 1] * z[ip1 - 1] + delta[i - 1] * delta[i - 1] * (dpsi + dphi);
                            }
                        } else {
                            a = delta[i - 1] * delta[i - 1] * dpsi + delta[ip1 - 1] * delta[ip1 - 1] * dphi;
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
                temp = rhoinv + psi + phi;
                if (swtch) {
                    c = temp - delta[iim1 - 1] * dpsi - delta[iip1 - 1] * dphi;
                    zz[1 - 1] = delta[iim1 - 1] * delta[iim1 - 1] * dpsi;
                    zz[3 - 1] = delta[iip1 - 1] * delta[iip1 - 1] * dphi;
                } else {
                    if (orgati) {
                        temp1 = z[iim1 - 1] / delta[iim1 - 1];
                        temp1 = temp1 * temp1;
                        c = temp - delta[iip1 - 1] * (dpsi + dphi) - (d[iim1 - 1] - d[iip1 - 1]) * temp1;
                        zz[1 - 1] = z[iim1 - 1] * z[iim1 - 1];
                        zz[3 - 1] = delta[iip1 - 1] * delta[iip1 - 1] * ((dpsi - temp1) + dphi);
                    } else {
                        temp1 = z[iip1 - 1] / delta[iip1 - 1];
                        temp1 = temp1 * temp1;
                        c = temp - delta[iim1 - 1] * (dpsi + dphi) - (d[iip1 - 1] - d[iim1 - 1]) * temp1;
                        zz[1 - 1] = delta[iim1 - 1] * delta[iim1 - 1] * (dpsi + (dphi - temp1));
                        zz[3 - 1] = z[iip1 - 1] * z[iip1 - 1];
                    }
                }
                Rlaed6(niter, orgati, c, delta[iim1 - 1], zz, w, eta, info);
                if (info != 0) {
                    goto statement_250;
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
            temp = tau + eta;
            if (temp > dltub || temp < dltlb) {
                if (w < zero) {
                    eta = (dltub - tau) / two;
                } else {
                    eta = (dltlb - tau) / two;
                }
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                delta[j - 1] = delta[j - 1] - eta;
            }
            //
            tau += eta;
            prew = w;
            //
            //           Evaluate PSI and the derivative DPSI
            //
            dpsi = zero;
            psi = zero;
            erretm = zero;
            for (j = 1; j <= iim1; j = j + 1) {
                temp = z[j - 1] / delta[j - 1];
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
                temp = z[j - 1] / delta[j - 1];
                phi += z[j - 1] * temp;
                dphi += temp * temp;
                erretm += phi;
            }
            //
            temp = z[ii - 1] / delta[ii - 1];
            dw = dpsi + dphi + temp * temp;
            temp = z[ii - 1] * temp;
            w = rhoinv + phi + psi + temp;
            erretm = eight * (phi - psi) + erretm + two * rhoinv + three * abs(temp) + abs(tau) * dw;
            if (w * prew > zero && abs(w) > abs(prew) / ten) {
                swtch = !swtch;
            }
            //
        }
        //
        //        Return with INFO = 1, NITER = MAXIT and not converged
        //
        info = 1;
        if (orgati) {
            dlam = d[i - 1] + tau;
        } else {
            dlam = d[ip1 - 1] + tau;
        }
        //
    }
//
statement_250:;
    //
    //     End of Rlaed4
    //
}
