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

void Rlahqr(bool const wantt, bool const wantz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, REAL *h, INTEGER const ldh, REAL *wr, REAL *wi, INTEGER const iloz, INTEGER const ihiz, REAL *z, INTEGER const ldz, INTEGER &info) {
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER nh = 0;
    INTEGER nz = 0;
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER itmax = 0;
    INTEGER kdefl = 0;
    INTEGER i = 0;
    INTEGER l = 0;
    INTEGER its = 0;
    INTEGER k = 0;
    REAL tst = 0.0;
    REAL ab = 0.0;
    REAL ba = 0.0;
    REAL aa = 0.0;
    REAL bb = 0.0;
    REAL s = 0.0;
    const INTEGER kexsh = 10;
    const REAL dat1 = 3.0 / 4.0;
    REAL h11 = 0.0;
    const REAL dat2 = -0.4375e0;
    REAL h12 = 0.0;
    REAL h21 = 0.0;
    REAL h22 = 0.0;
    REAL rt1r = 0.0;
    REAL rt1i = 0.0;
    REAL rt2r = 0.0;
    REAL rt2i = 0.0;
    const REAL two = 2.0;
    REAL tr = 0.0;
    REAL det = 0.0;
    REAL rtdisc = 0.0;
    INTEGER m = 0;
    REAL h21s = 0.0;
    REAL v[3];
    INTEGER nr = 0;
    REAL t1 = 0.0;
    REAL v2 = 0.0;
    REAL t2 = 0.0;
    REAL v3 = 0.0;
    REAL t3 = 0.0;
    REAL sum = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    if (ilo == ihi) {
        wr[ilo - 1] = h[(ilo - 1) + (ilo - 1) * ldh];
        wi[ilo - 1] = zero;
        return;
    }
    //
    //     ==== clear out the trash ====
    for (j = ilo; j <= ihi - 3; j = j + 1) {
        h[((j + 2) - 1) + (j - 1) * ldh] = zero;
        h[((j + 3) - 1) + (j - 1) * ldh] = zero;
    }
    if (ilo <= ihi - 2) {
        h[(ihi - 1) + ((ihi - 2) - 1) * ldh] = zero;
    }
    //
    nh = ihi - ilo + 1;
    nz = ihiz - iloz + 1;
    //
    //     Set machine-dependent constants for the stopping criterion.
    //
    safmin = Rlamch("SAFE MINIMUM");
    safmax = one / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * (castREAL(nh) / ulp);
    //
    //     I1 and I2 are the indices of the first row and last column of H
    //     to which transformations must be applied. If eigenvalues only are
    //     being computed, I1 and I2 are set inside the main loop.
    //
    if (wantt) {
        i1 = 1;
        i2 = n;
    }
    //
    //     ITMAX is the total number of QR iterations allowed.
    //
    itmax = 30 * max((INTEGER)10, nh);
    //
    //     KDEFL counts the number of iterations since a deflation
    //
    kdefl = 0;
    //
    //     The main loop begins here. I is the loop index and decreases from
    //     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
    //     with the active submatrix in rows and columns L to I.
    //     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
    //     H(L,L-1) is negligible so that the matrix splits.
    //
    i = ihi;
statement_20:
    l = ilo;
    if (i < ilo) {
        goto statement_160;
    }
    //
    //     Perform QR iterations on rows and columns ILO to I until a
    //     submatrix of order 1 or 2 splits off at the bottom because a
    //     subdiagonal element has become negligible.
    //
    for (its = 0; its <= itmax; its = its + 1) {
        //
        //        Look for a single small subdiagonal element.
        //
        for (k = i; k >= l + 1; k = k - 1) {
            if (abs(h[(k - 1) + ((k - 1) - 1) * ldh]) <= smlnum) {
                goto statement_40;
            }
            tst = abs(h[((k - 1) - 1) + ((k - 1) - 1) * ldh]) + abs(h[(k - 1) + (k - 1) * ldh]);
            if (tst == zero) {
                if (k - 2 >= ilo) {
                    tst += abs(h[((k - 1) - 1) + ((k - 2) - 1) * ldh]);
                }
                if (k + 1 <= ihi) {
                    tst += abs(h[((k + 1) - 1) + (k - 1) * ldh]);
                }
            }
            //           ==== The following is a conservative small subdiagonal
            //           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
            //           .    1997). It has better mathematical foundation and
            //           .    improves accuracy in some cases.  ====
            if (abs(h[(k - 1) + ((k - 1) - 1) * ldh]) <= ulp * tst) {
                ab = max(abs(h[(k - 1) + ((k - 1) - 1) * ldh]), abs(h[((k - 1) - 1) + (k - 1) * ldh]));
                ba = min(abs(h[(k - 1) + ((k - 1) - 1) * ldh]), abs(h[((k - 1) - 1) + (k - 1) * ldh]));
                aa = max(REAL(abs(h[(k - 1) + (k - 1) * ldh])), REAL(abs(h[((k - 1) - 1) + ((k - 1) - 1) * ldh] - h[(k - 1) + (k - 1) * ldh])));
                bb = min(REAL(abs(h[(k - 1) + (k - 1) * ldh])), REAL(abs(h[((k - 1) - 1) + ((k - 1) - 1) * ldh] - h[(k - 1) + (k - 1) * ldh])));
                s = aa + ab;
                if (ba * (ab / s) <= max(smlnum, REAL(ulp * (bb * (aa / s))))) {
                    goto statement_40;
                }
            }
        }
    statement_40:
        l = k;
        if (l > ilo) {
            //
            //           H(L,L-1) is negligible
            //
            h[(l - 1) + ((l - 1) - 1) * ldh] = zero;
        }
        //
        //        Exit from loop if a submatrix of order 1 or 2 has split off.
        //
        if (l >= i - 1) {
            goto statement_150;
        }
        kdefl++;
        //
        //        Now the active submatrix is in rows and columns L to I. If
        //        eigenvalues only are being computed, only the active submatrix
        //        need be transformed.
        //
        if (!wantt) {
            i1 = l;
            i2 = i;
        }
        //
        if (mod(kdefl, 2 * kexsh) == 0) {
            //
            //           Exceptional shift.
            //
            s = abs(h[(i - 1) + ((i - 1) - 1) * ldh]) + abs(h[((i - 1) - 1) + ((i - 2) - 1) * ldh]);
            h11 = dat1 * s + h[(i - 1) + (i - 1) * ldh];
            h12 = dat2 * s;
            h21 = s;
            h22 = h11;
        } else if (mod(kdefl, kexsh) == 0) {
            //
            //           Exceptional shift.
            //
            s = abs(h[((l + 1) - 1) + (l - 1) * ldh]) + abs(h[((l + 2) - 1) + ((l + 1) - 1) * ldh]);
            h11 = dat1 * s + h[(l - 1) + (l - 1) * ldh];
            h12 = dat2 * s;
            h21 = s;
            h22 = h11;
        } else {
            //
            //           Prepare to use Francis' REAL shift
            //           (i.e. 2nd degree generalized Rayleigh quotient)
            //
            h11 = h[((i - 1) - 1) + ((i - 1) - 1) * ldh];
            h21 = h[(i - 1) + ((i - 1) - 1) * ldh];
            h12 = h[((i - 1) - 1) + (i - 1) * ldh];
            h22 = h[(i - 1) + (i - 1) * ldh];
        }
        s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
        if (s == zero) {
            rt1r = zero;
            rt1i = zero;
            rt2r = zero;
            rt2i = zero;
        } else {
            h11 = h11 / s;
            h21 = h21 / s;
            h12 = h12 / s;
            h22 = h22 / s;
            tr = (h11 + h22) / two;
            det = (h11 - tr) * (h22 - tr) - h12 * h21;
            rtdisc = sqrt(abs(det));
            if (det >= zero) {
                //
                //              ==== complex conjugate shifts ====
                //
                rt1r = tr * s;
                rt2r = rt1r;
                rt1i = rtdisc * s;
                rt2i = -rt1i;
            } else {
                //
                //              ==== real shifts (use only one of them)  ====
                //
                rt1r = tr + rtdisc;
                rt2r = tr - rtdisc;
                if (abs(rt1r - h22) <= abs(rt2r - h22)) {
                    rt1r = rt1r * s;
                    rt2r = rt1r;
                } else {
                    rt2r = rt2r * s;
                    rt1r = rt2r;
                }
                rt1i = zero;
                rt2i = zero;
            }
        }
        //
        //        Look for two consecutive small subdiagonal elements.
        //
        for (m = i - 2; m >= l; m = m - 1) {
            //           Determine the effect of starting the REAL-shift QR
            //           iteration at row M, and see if this would make H(M,M-1)
            //           negligible.  (The following uses scaling to avoid
            //           overflows and most underflows.)
            //
            h21s = h[((m + 1) - 1) + (m - 1) * ldh];
            s = abs(h[(m - 1) + (m - 1) * ldh] - rt2r) + abs(rt2i) + abs(h21s);
            h21s = h[((m + 1) - 1) + (m - 1) * ldh] / s;
            v[1 - 1] = h21s * h[(m - 1) + ((m + 1) - 1) * ldh] + (h[(m - 1) + (m - 1) * ldh] - rt1r) * ((h[(m - 1) + (m - 1) * ldh] - rt2r) / s) - rt1i * (rt2i / s);
            v[2 - 1] = h21s * (h[(m - 1) + (m - 1) * ldh] + h[((m + 1) - 1) + ((m + 1) - 1) * ldh] - rt1r - rt2r);
            v[3 - 1] = h21s * h[((m + 2) - 1) + ((m + 1) - 1) * ldh];
            s = abs(v[1 - 1]) + abs(v[2 - 1]) + abs(v[3 - 1]);
            v[1 - 1] = v[1 - 1] / s;
            v[2 - 1] = v[2 - 1] / s;
            v[3 - 1] = v[3 - 1] / s;
            if (m == l) {
                goto statement_60;
            }
            if (abs(h[(m - 1) + ((m - 1) - 1) * ldh]) * (abs(v[2 - 1]) + abs(v[3 - 1])) <= ulp * abs(v[1 - 1]) * (abs(h[((m - 1) - 1) + ((m - 1) - 1) * ldh]) + abs(h[(m - 1) + (m - 1) * ldh]) + abs(h[((m + 1) - 1) + ((m + 1) - 1) * ldh]))) {
                goto statement_60;
            }
        }
    statement_60:
        //
        //        Double-shift QR step
        //
        for (k = m; k <= i - 1; k = k + 1) {
            //
            //           The first iteration of this loop determines a reflection G
            //           from the vector V and applies it from left and right to H,
            //           thus creating a nonzero bulge below the subdiagonal.
            //
            //           Each subsequent iteration determines a reflection G to
            //           restore the Hessenberg form in the (K-1)th column, and thus
            //           chases the bulge one step toward the bottom of the active
            //           submatrix. NR is the order of G.
            //
            nr = min((INTEGER)3, i - k + 1);
            if (k > m) {
                Rcopy(nr, &h[(k - 1) + ((k - 1) - 1) * ldh], 1, v, 1);
            }
            Rlarfg(nr, v[1 - 1], &v[2 - 1], 1, t1);
            if (k > m) {
                h[(k - 1) + ((k - 1) - 1) * ldh] = v[1 - 1];
                h[((k + 1) - 1) + ((k - 1) - 1) * ldh] = zero;
                if (k < i - 1) {
                    h[((k + 2) - 1) + ((k - 1) - 1) * ldh] = zero;
                }
            } else if (m > l) {
                //               ==== Use the following instead of
                //               .    H( K, K-1 ) = -H( K, K-1 ) to
                //               .    avoid a bug when v(2) and v(3)
                //               .    underflow. ====
                h[(k - 1) + ((k - 1) - 1) * ldh] = h[(k - 1) + ((k - 1) - 1) * ldh] * (one - t1);
            }
            v2 = v[2 - 1];
            t2 = t1 * v2;
            if (nr == 3) {
                v3 = v[3 - 1];
                t3 = t1 * v3;
                //
                //              Apply G from the left to transform the rows of the matrix
                //              in columns K to I2.
                //
                for (j = k; j <= i2; j = j + 1) {
                    sum = h[(k - 1) + (j - 1) * ldh] + v2 * h[((k + 1) - 1) + (j - 1) * ldh] + v3 * h[((k + 2) - 1) + (j - 1) * ldh];
                    h[(k - 1) + (j - 1) * ldh] = h[(k - 1) + (j - 1) * ldh] - sum * t1;
                    h[((k + 1) - 1) + (j - 1) * ldh] = h[((k + 1) - 1) + (j - 1) * ldh] - sum * t2;
                    h[((k + 2) - 1) + (j - 1) * ldh] = h[((k + 2) - 1) + (j - 1) * ldh] - sum * t3;
                }
                //
                //              Apply G from the right to transform the columns of the
                //              matrix in rows I1 to min(K+3,I).
                //
                for (j = i1; j <= min(k + 3, i); j = j + 1) {
                    sum = h[(j - 1) + (k - 1) * ldh] + v2 * h[(j - 1) + ((k + 1) - 1) * ldh] + v3 * h[(j - 1) + ((k + 2) - 1) * ldh];
                    h[(j - 1) + (k - 1) * ldh] = h[(j - 1) + (k - 1) * ldh] - sum * t1;
                    h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - sum * t2;
                    h[(j - 1) + ((k + 2) - 1) * ldh] = h[(j - 1) + ((k + 2) - 1) * ldh] - sum * t3;
                }
                //
                if (wantz) {
                    //
                    //                 Accumulate transformations in the matrix Z
                    //
                    for (j = iloz; j <= ihiz; j = j + 1) {
                        sum = z[(j - 1) + (k - 1) * ldz] + v2 * z[(j - 1) + ((k + 1) - 1) * ldz] + v3 * z[(j - 1) + ((k + 2) - 1) * ldz];
                        z[(j - 1) + (k - 1) * ldz] = z[(j - 1) + (k - 1) * ldz] - sum * t1;
                        z[(j - 1) + ((k + 1) - 1) * ldz] = z[(j - 1) + ((k + 1) - 1) * ldz] - sum * t2;
                        z[(j - 1) + ((k + 2) - 1) * ldz] = z[(j - 1) + ((k + 2) - 1) * ldz] - sum * t3;
                    }
                }
            } else if (nr == 2) {
                //
                //              Apply G from the left to transform the rows of the matrix
                //              in columns K to I2.
                //
                for (j = k; j <= i2; j = j + 1) {
                    sum = h[(k - 1) + (j - 1) * ldh] + v2 * h[((k + 1) - 1) + (j - 1) * ldh];
                    h[(k - 1) + (j - 1) * ldh] = h[(k - 1) + (j - 1) * ldh] - sum * t1;
                    h[((k + 1) - 1) + (j - 1) * ldh] = h[((k + 1) - 1) + (j - 1) * ldh] - sum * t2;
                }
                //
                //              Apply G from the right to transform the columns of the
                //              matrix in rows I1 to min(K+3,I).
                //
                for (j = i1; j <= i; j = j + 1) {
                    sum = h[(j - 1) + (k - 1) * ldh] + v2 * h[(j - 1) + ((k + 1) - 1) * ldh];
                    h[(j - 1) + (k - 1) * ldh] = h[(j - 1) + (k - 1) * ldh] - sum * t1;
                    h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - sum * t2;
                }
                //
                if (wantz) {
                    //
                    //                 Accumulate transformations in the matrix Z
                    //
                    for (j = iloz; j <= ihiz; j = j + 1) {
                        sum = z[(j - 1) + (k - 1) * ldz] + v2 * z[(j - 1) + ((k + 1) - 1) * ldz];
                        z[(j - 1) + (k - 1) * ldz] = z[(j - 1) + (k - 1) * ldz] - sum * t1;
                        z[(j - 1) + ((k + 1) - 1) * ldz] = z[(j - 1) + ((k + 1) - 1) * ldz] - sum * t2;
                    }
                }
            }
        }
        //
    }
    //
    //     Failure to converge in remaining number of iterations
    //
    info = i;
    return;
//
statement_150:
    //
    if (l == i) {
        //
        //        H(I,I-1) is negligible: one eigenvalue has converged.
        //
        wr[i - 1] = h[(i - 1) + (i - 1) * ldh];
        wi[i - 1] = zero;
    } else if (l == i - 1) {
        //
        //        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
        //
        //        Transform the 2-by-2 submatrix to standard Schur form,
        //        and compute and store the eigenvalues.
        //
        Rlanv2(h[((i - 1) - 1) + ((i - 1) - 1) * ldh], h[((i - 1) - 1) + (i - 1) * ldh], h[(i - 1) + ((i - 1) - 1) * ldh], h[(i - 1) + (i - 1) * ldh], wr[(i - 1) - 1], wi[(i - 1) - 1], wr[i - 1], wi[i - 1], cs, sn);
        //
        if (wantt) {
            //
            //           Apply the transformation to the rest of H.
            //
            if (i2 > i) {
                Rrot(i2 - i, &h[((i - 1) - 1) + ((i + 1) - 1) * ldh], ldh, &h[(i - 1) + ((i + 1) - 1) * ldh], ldh, cs, sn);
            }
            Rrot(i - i1 - 1, &h[(i1 - 1) + ((i - 1) - 1) * ldh], 1, &h[(i1 - 1) + (i - 1) * ldh], 1, cs, sn);
        }
        if (wantz) {
            //
            //           Apply the transformation to Z.
            //
            Rrot(nz, &z[(iloz - 1) + ((i - 1) - 1) * ldz], 1, &z[(iloz - 1) + (i - 1) * ldz], 1, cs, sn);
        }
    }
    //     reset deflation counter
    kdefl = 0;
    //
    //     return to start of the main loop with new value of I.
    //
    i = l - 1;
    goto statement_20;
//
statement_160:;
    //
    //     End of Rlahqr
    //
}
