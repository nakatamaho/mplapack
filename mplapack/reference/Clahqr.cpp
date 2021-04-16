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

inline REAL abs1(COMPLEX ff) { return max(abs(ff.real()), abs(ff.imag())); }

void Clahqr(bool const wantt, bool const wantz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, COMPLEX *h, INTEGER const ldh, COMPLEX *w, INTEGER const iloz, INTEGER const ihiz, COMPLEX *z, INTEGER const ldz, INTEGER &info) {
    COMPLEX cdum = 0.0;
    INTEGER j = 0;
    const COMPLEX zero = (0.0, 0.0);
    INTEGER jlo = 0;
    INTEGER jhi = 0;
    INTEGER i = 0;
    const REAL rzero = 0.0;
    COMPLEX sc = 0.0;
    INTEGER nh = 0;
    INTEGER nz = 0;
    REAL safmin = 0.0;
    const REAL rone = 1.0;
    REAL safmax = 0.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER itmax = 0;
    INTEGER kdefl = 0;
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
    COMPLEX t = 0.0;
    COMPLEX u = 0.0;
    const REAL half = 0.5e0;
    COMPLEX x = 0.0;
    REAL sx = 0.0;
    COMPLEX y = 0.0;
    INTEGER m = 0;
    COMPLEX h11 = 0.0;
    COMPLEX h22 = 0.0;
    COMPLEX h11s = 0.0;
    REAL h21 = 0.0;
    COMPLEX v[2];
    REAL h10 = 0.0;
    COMPLEX t1 = 0.0;
    COMPLEX v2 = 0.0;
    REAL t2 = 0.0;
    COMPLEX sum = 0.0;
    const COMPLEX one = (1.0, 0.0);
    COMPLEX temp = 0.0;
    REAL rtemp = 0.0;
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
    //  =========================================================
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
    //     .. Statement Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Function definitions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    if (ilo == ihi) {
        w[ilo - 1] = h[(ilo - 1) + (ilo - 1) * ldh];
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
    //     ==== ensure that subdiagonal entries are real ====
    if (wantt) {
        jlo = 1;
        jhi = n;
    } else {
        jlo = ilo;
        jhi = ihi;
    }
    for (i = ilo + 1; i <= ihi; i = i + 1) {
        if (h[(i - 1) + ((i - 1) - 1) * ldh].imag() != rzero) {
            //           ==== The following redundant normalization
            //           .    avoids problems with both gradual and
            //           .    sudden underflow in ABS(H(I,I-1)) ====
            sc = h[(i - 1) + ((i - 1) - 1) * ldh] / abs1(h[(i - 1) + ((i - 1) - 1) * ldh]);
            sc = conj(sc) / abs(sc);
            h[(i - 1) + ((i - 1) - 1) * ldh] = abs(h[(i - 1) + ((i - 1) - 1) * ldh]);
            Cscal(jhi - i + 1, sc, &h[(i - 1) + (i - 1) * ldh], ldh);
            Cscal(min(jhi, i + 1) - jlo + 1, conj(sc), &h[(jlo - 1) + (i - 1) * ldh], 1);
            if (wantz) {
                Cscal(ihiz - iloz + 1, conj(sc), &z[(iloz - 1) + (i - 1) * ldz], 1);
            }
        }
    }
    //
    nh = ihi - ilo + 1;
    nz = ihiz - iloz + 1;
    //
    //     Set machine-dependent constants for the stopping criterion.
    //
    safmin = Rlamch("SAFE MINIMUM");
    safmax = rone / safmin;
    Rlabad(safmin, safmax);
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
    //     IHI to ILO in steps of 1. Each iteration of the loop works
    //     with the active submatrix in rows and columns L to I.
    //     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
    //     H(L,L-1) is negligible so that the matrix splits.
    //
    i = ihi;
statement_30:
    if (i < ilo) {
        goto statement_150;
    }
    //
    //     Perform QR iterations on rows and columns ILO to I until a
    //     submatrix of order 1 splits off at the bottom because a
    //     subdiagonal element has become negligible.
    //
    l = ilo;
    for (its = 0; its <= itmax; its = its + 1) {
        //
        //        Look for a single small subdiagonal element.
        //
        for (k = i; k >= l + 1; k = k - 1) {
            if (abs1(h[(k - 1) + ((k - 1) - 1) * ldh]) <= smlnum) {
                goto statement_50;
            }
            tst = abs1(h[((k - 1) - 1) + ((k - 1) - 1) * ldh]) + abs1(h[(k - 1) + (k - 1) * ldh]);
            if (tst == zero) {
                if (k - 2 >= ilo) {
                    tst += abs(h[((k - 1) - 1) + ((k - 2) - 1) * ldh].real());
                }
                if (k + 1 <= ihi) {
                    tst += abs(h[((k + 1) - 1) + (k - 1) * ldh].real());
                }
            }
            //           ==== The following is a conservative small subdiagonal
            //           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
            //           .    1997). It has better mathematical foundation and
            //           .    improves accuracy in some examples.  ====
            if (abs(h[(k - 1) + ((k - 1) - 1) * ldh].real()) <= ulp * tst) {
                ab = max(abs1(h[(k - 1) + ((k - 1) - 1) * ldh]), abs1(h[((k - 1) - 1) + (k - 1) * ldh]));
                ba = min(abs1(h[(k - 1) + ((k - 1) - 1) * ldh]), abs1(h[((k - 1) - 1) + (k - 1) * ldh]));
                aa = max(abs1(h[(k - 1) + (k - 1) * ldh]), abs1(h[((k - 1) - 1) + ((k - 1) - 1) * ldh] - h[(k - 1) + (k - 1) * ldh]));
                bb = min(abs1(h[(k - 1) + (k - 1) * ldh]), abs1(h[((k - 1) - 1) + ((k - 1) - 1) * ldh] - h[(k - 1) + (k - 1) * ldh]));
                s = aa + ab;
                if (ba * (ab / s) <= max(smlnum, ulp * (bb * (aa / s)))) {
                    goto statement_50;
                }
            }
        }
    statement_50:
        l = k;
        if (l > ilo) {
            //
            //           H(L,L-1) is negligible
            //
            h[(l - 1) + ((l - 1) - 1) * ldh] = zero;
        }
        //
        //        Exit from loop if a submatrix of order 1 has split off.
        //
        if (l >= i) {
            goto statement_140;
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
            s = dat1 * abs(h[(i - 1) + ((i - 1) - 1) * ldh].real());
            t = s + h[(i - 1) + (i - 1) * ldh];
        } else if (mod(kdefl, kexsh) == 0) {
            //
            //           Exceptional shift.
            //
            s = dat1 * abs(h[((l + 1) - 1) + (l - 1) * ldh].real());
            t = s + h[(l - 1) + (l - 1) * ldh];
        } else {
            //
            //           Wilkinson's shift.
            //
            t = h[(i - 1) + (i - 1) * ldh];
            u = sqrt(h[((i - 1) - 1) + (i - 1) * ldh]) * sqrt(h[(i - 1) + ((i - 1) - 1) * ldh]);
            s = abs1(u);
            if (s != rzero) {
                x = half * (h[((i - 1) - 1) + ((i - 1) - 1) * ldh] - t);
                sx = abs1(x);
                s = max(s, abs1(x));
                y = s * sqrt((x / s) * (x / s) + (u / s) * (u / s));
                if (sx > rzero) {
                    if ((x / sx).real() * y.real() + (x / sx).imag() * y.imag() < rzero) {
                        y = -y;
                    }
                }
                t = t - u * Cladiv(u, (x + y));
            }
        }
        //
        //        Look for two consecutive small subdiagonal elements.
        //
        for (m = i - 1; m >= l + 1; m = m - 1) {
            //
            //           Determine the effect of starting the single-shift QR
            //           iteration at row M, and see if this would make H(M,M-1)
            //           negligible.
            //
            h11 = h[(m - 1) + (m - 1) * ldh];
            h22 = h[((m + 1) - 1) + ((m + 1) - 1) * ldh];
            h11s = h11 - t;
            h21 = h[((m + 1) - 1) + (m - 1) * ldh].real();
            s = abs1(h11s) + abs(h21);
            h11s = h11s / s;
            h21 = h21 / s;
            v[1 - 1] = h11s;
            v[2 - 1] = h21;
            h10 = h[(m - 1) + ((m - 1) - 1) * ldh].real();
            if (abs(h10) * abs(h21) <= ulp * (abs1(h11s) * (abs1(h11) + abs1(h22)))) {
                goto statement_70;
            }
        }
        h11 = h[(l - 1) + (l - 1) * ldh];
        h22 = h[((l + 1) - 1) + ((l + 1) - 1) * ldh];
        h11s = h11 - t;
        h21 = h[((l + 1) - 1) + (l - 1) * ldh].real();
        s = abs1(h11s) + abs(h21);
        h11s = h11s / s;
        h21 = h21 / s;
        v[1 - 1] = h11s;
        v[2 - 1] = h21;
    statement_70:
        //
        //        Single-shift QR step
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
            //           submatrix.
            //
            //           V(2) is always real before the call to Clarfg, and hence
            //           after the call T2 ( = T1*V(2) ) is also real.
            //
            if (k > m) {
                Ccopy(2, &h[(k - 1) + ((k - 1) - 1) * ldh], 1, v, 1);
            }
            Clarfg(2, v[1 - 1], &v[2 - 1], 1, t1);
            if (k > m) {
                h[(k - 1) + ((k - 1) - 1) * ldh] = v[1 - 1];
                h[((k + 1) - 1) + ((k - 1) - 1) * ldh] = zero;
            }
            v2 = v[2 - 1];
            t2 = (t1 * v2).real();
            //
            //           Apply G from the left to transform the rows of the matrix
            //           in columns K to I2.
            //
            for (j = k; j <= i2; j = j + 1) {
                sum = conj(t1) * h[(k - 1) + (j - 1) * ldh] + t2 * h[((k + 1) - 1) + (j - 1) * ldh];
                h[(k - 1) + (j - 1) * ldh] = h[(k - 1) + (j - 1) * ldh] - sum;
                h[((k + 1) - 1) + (j - 1) * ldh] = h[((k + 1) - 1) + (j - 1) * ldh] - sum * v2;
            }
            //
            //           Apply G from the right to transform the columns of the
            //           matrix in rows I1 to min(K+2,I).
            //
            for (j = i1; j <= min(k + 2, i); j = j + 1) {
                sum = t1 * h[(j - 1) + (k - 1) * ldh] + t2 * h[(j - 1) + ((k + 1) - 1) * ldh];
                h[(j - 1) + (k - 1) * ldh] = h[(j - 1) + (k - 1) * ldh] - sum;
                h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - sum * conj(v2);
            }
            //
            if (wantz) {
                //
                //              Accumulate transformations in the matrix Z
                //
                for (j = iloz; j <= ihiz; j = j + 1) {
                    sum = t1 * z[(j - 1) + (k - 1) * ldz] + t2 * z[(j - 1) + ((k + 1) - 1) * ldz];
                    z[(j - 1) + (k - 1) * ldz] = z[(j - 1) + (k - 1) * ldz] - sum;
                    z[(j - 1) + ((k + 1) - 1) * ldz] = z[(j - 1) + ((k + 1) - 1) * ldz] - sum * conj(v2);
                }
            }
            //
            if (k == m && m > l) {
                //
                //              If the QR step was started at row M > L because two
                //              consecutive small subdiagonals were found, then extra
                //              scaling must be performed to ensure that H(M,M-1) remains
                //              real.
                //
                temp = one - t1;
                temp = temp / abs(temp);
                h[((m + 1) - 1) + (m - 1) * ldh] = h[((m + 1) - 1) + (m - 1) * ldh] * conj(temp);
                if (m + 2 <= i) {
                    h[((m + 2) - 1) + ((m + 1) - 1) * ldh] = h[((m + 2) - 1) + ((m + 1) - 1) * ldh] * temp;
                }
                for (j = m; j <= i; j = j + 1) {
                    if (j != m + 1) {
                        if (i2 > j) {
                            Cscal(i2 - j, temp, &h[(j - 1) + ((j + 1) - 1) * ldh], ldh);
                        }
                        Cscal(j - i1, conj(temp), &h[(i1 - 1) + (j - 1) * ldh], 1);
                        if (wantz) {
                            Cscal(nz, conj(temp), &z[(iloz - 1) + (j - 1) * ldz], 1);
                        }
                    }
                }
            }
        }
        //
        //        Ensure that H(I,I-1) is real.
        //
        temp = h[(i - 1) + ((i - 1) - 1) * ldh];
        if (temp.imag() != rzero) {
            rtemp = abs(temp);
            h[(i - 1) + ((i - 1) - 1) * ldh] = rtemp;
            temp = temp / rtemp;
            if (i2 > i) {
                Cscal(i2 - i, conj(temp), &h[(i - 1) + ((i + 1) - 1) * ldh], ldh);
            }
            Cscal(i - i1, temp, &h[(i1 - 1) + (i - 1) * ldh], 1);
            if (wantz) {
                Cscal(nz, temp, &z[(iloz - 1) + (i - 1) * ldz], 1);
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
statement_140:
    //
    //     H(I,I-1) is negligible: one eigenvalue has converged.
    //
    w[i - 1] = h[(i - 1) + (i - 1) * ldh];
    //     reset deflation counter
    kdefl = 0;
    //
    //     return to start of the main loop with new value of I.
    //
    i = l - 1;
    goto statement_30;
//
statement_150:;
    //
    //     End of Clahqr
    //
}
