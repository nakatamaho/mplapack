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

void Rlasq2(INTEGER const n, REAL *z, INTEGER &info) {
    REAL eps = 0.0;
    REAL safmin = 0.0;
    const REAL hundrd = 100.0;
    REAL tol = 0.0;
    REAL tol2 = 0.0;
    const REAL zero = 0.0;
    REAL d = 0.0;
    const REAL half = 0.5e0;
    REAL t = 0.0;
    REAL s = 0.0;
    const REAL one = 1.0;
    REAL emin = 0.0;
    REAL qmax = 0.0;
    REAL zmax = 0.0;
    REAL e = 0.0;
    INTEGER k = 0;
    INTEGER iinfo = 0;
    REAL trace = 0.0;
    bool ieee = false;
    INTEGER i0 = 0;
    INTEGER n0 = 0;
    const REAL cbias = 1.50e0;
    INTEGER ipn4 = 0;
    INTEGER i4 = 0;
    REAL temp = 0.0;
    INTEGER pp = 0;
    INTEGER ttype = 0;
    REAL dmin1 = 0.0;
    REAL dmin2 = 0.0;
    REAL dn = 0.0;
    REAL dn1 = 0.0;
    REAL dn2 = 0.0;
    REAL g = 0.0;
    REAL tau = 0.0;
    INTEGER iter = 0;
    INTEGER nfail = 0;
    INTEGER ndiv = 0;
    INTEGER iwhila = 0;
    REAL desig = 0.0;
    REAL sigma = 0.0;
    REAL emax = 0.0;
    REAL qmin = 0.0;
    const REAL four = 4.0;
    REAL dee = 0.0;
    REAL deemin = 0.0;
    INTEGER kmin = 0;
    const REAL two = 2.0;
    REAL dmin = 0.0;
    INTEGER nbig = 0;
    INTEGER iwhilb = 0;
    INTEGER splt = 0;
    REAL oldemn = 0.0;
    INTEGER i1 = 0;
    INTEGER n1 = 0;
    REAL tempq = 0.0;
    REAL tempe = 0.0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments.
    //     (in case Rlasq2 is not called by Rlasq1)
    //
    info = 0;
    eps = Rlamch("Precision");
    safmin = Rlamch("Safe minimum");
    tol = eps * hundrd;
    tol2 = pow2(tol);
    //
    if (n < 0) {
        info = -1;
        Mxerbla("Rlasq2", 1);
        return;
    } else if (n == 0) {
        return;
    } else if (n == 1) {
        //
        //        1-by-1 case.
        //
        if (z[1 - 1] < zero) {
            info = -201;
            Mxerbla("Rlasq2", 2);
        }
        return;
    } else if (n == 2) {
        //
        //        2-by-2 case.
        //
        if (z[1 - 1] < zero) {
            info = -201;
            Mxerbla("Rlasq2", 2);
            return;
        } else if (z[2 - 1] < zero) {
            info = -202;
            Mxerbla("Rlasq2", 2);
            return;
        } else if (z[3 - 1] < zero) {
            info = -203;
            Mxerbla("Rlasq2", 2);
            return;
        } else if (z[3 - 1] > z[1 - 1]) {
            d = z[3 - 1];
            z[3 - 1] = z[1 - 1];
            z[1 - 1] = d;
        }
        z[5 - 1] = z[1 - 1] + z[2 - 1] + z[3 - 1];
        if (z[2 - 1] > z[3 - 1] * tol2) {
            t = half * ((z[1 - 1] - z[3 - 1]) + z[2 - 1]);
            s = z[3 - 1] * (z[2 - 1] / t);
            if (s <= t) {
                s = z[3 - 1] * (z[2 - 1] / (t * (one + sqrt(one + s / t))));
            } else {
                s = z[3 - 1] * (z[2 - 1] / (t + sqrt(t) * sqrt(t + s)));
            }
            t = z[1 - 1] + (s + z[2 - 1]);
            z[3 - 1] = z[3 - 1] * (z[1 - 1] / t);
            z[1 - 1] = t;
        }
        z[2 - 1] = z[3 - 1];
        z[6 - 1] = z[2 - 1] + z[1 - 1];
        return;
    }
    //
    //     Check for negative data and compute sums of q's and e's.
    //
    z[(2 * n) - 1] = zero;
    emin = z[2 - 1];
    qmax = zero;
    zmax = zero;
    d = zero;
    e = zero;
    //
    for (k = 1; k <= 2 * (n - 1); k = k + 2) {
        if (z[k - 1] < zero) {
            info = -(200 + k);
            Mxerbla("Rlasq2", 2);
            return;
        } else if (z[(k + 1) - 1] < zero) {
            info = -(200 + k + 1);
            Mxerbla("Rlasq2", 2);
            return;
        }
        d += z[k - 1];
        e += z[(k + 1) - 1];
        qmax = max(qmax, &z[k - 1]);
        emin = min(emin, &z[(k + 1) - 1]);
        zmax = max(qmax, zmax, &z[(k + 1) - 1]);
    }
    if (z[(2 * n - 1) - 1] < zero) {
        info = -(200 + 2 * n - 1);
        Mxerbla("Rlasq2", 2);
        return;
    }
    d += z[(2 * n - 1) - 1];
    qmax = max(qmax, &z[(2 * n - 1) - 1]);
    zmax = max(qmax, zmax);
    //
    //     Check for diagonality.
    //
    if (e == zero) {
        for (k = 2; k <= n; k = k + 1) {
            z[k - 1] = z[(2 * k - 1) - 1];
        }
        Rlasrt("D", n, z, iinfo);
        z[(2 * n - 1) - 1] = d;
        return;
    }
    //
    trace = d + e;
    //
    //     Check for zero data.
    //
    if (trace == zero) {
        z[(2 * n - 1) - 1] = zero;
        return;
    }
    //
    //     Check whether the machine is IEEE conformable.
    //
    ieee = (iMlaenv(10, "Rlasq2", "N", 1, 2, 3, 4) == 1);
    //
    //     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
    //
    for (k = 2 * n; k >= 2; k = k - 2) {
        z[(2 * k) - 1] = zero;
        z[(2 * k - 1) - 1] = z[k - 1];
        z[(2 * k - 2) - 1] = zero;
        z[(2 * k - 3) - 1] = z[(k - 1) - 1];
    }
    //
    i0 = 1;
    n0 = n;
    //
    //     Reverse the qd-array, if warranted.
    //
    if (cbias * z[(4 * i0 - 3) - 1] < z[(4 * n0 - 3) - 1]) {
        ipn4 = 4 * (i0 + n0);
        for (i4 = 4 * i0; i4 <= 2 * (i0 + n0 - 1); i4 = i4 + 4) {
            temp = z[(i4 - 3) - 1];
            z[(i4 - 3) - 1] = z[(ipn4 - i4 - 3) - 1];
            z[(ipn4 - i4 - 3) - 1] = temp;
            temp = z[(i4 - 1) - 1];
            z[(i4 - 1) - 1] = z[(ipn4 - i4 - 5) - 1];
            z[(ipn4 - i4 - 5) - 1] = temp;
        }
    }
    //
    //     Initial split checking via dqd and Li's test.
    //
    pp = 0;
    //
    for (k = 1; k <= 2; k = k + 1) {
        //
        d = z[(4 * n0 + pp - 3) - 1];
        for (i4 = 4 * (n0 - 1) + pp; i4 >= 4 * i0 + pp; i4 = i4 - 4) {
            if (z[(i4 - 1) - 1] <= tol2 * d) {
                z[(i4 - 1) - 1] = -zero;
                d = z[(i4 - 3) - 1];
            } else {
                d = z[(i4 - 3) - 1] * (d / (d + z[(i4 - 1) - 1]));
            }
        }
        //
        //        dqd maps Z to ZZ plus Li's test.
        //
        emin = z[(4 * i0 + pp + 1) - 1];
        d = z[(4 * i0 + pp - 3) - 1];
        for (i4 = 4 * i0 + pp; i4 <= 4 * (n0 - 1) + pp; i4 = i4 + 4) {
            z[(i4 - 2 * pp - 2) - 1] = d + z[(i4 - 1) - 1];
            if (z[(i4 - 1) - 1] <= tol2 * d) {
                z[(i4 - 1) - 1] = -zero;
                z[(i4 - 2 * pp - 2) - 1] = d;
                z[(i4 - 2 * pp) - 1] = zero;
                d = z[(i4 + 1) - 1];
            } else if (safmin * z[(i4 + 1) - 1] < z[(i4 - 2 * pp - 2) - 1] && safmin * z[(i4 - 2 * pp - 2) - 1] < z[(i4 + 1) - 1]) {
                temp = z[(i4 + 1) - 1] / z[(i4 - 2 * pp - 2) - 1];
                z[(i4 - 2 * pp) - 1] = z[(i4 - 1) - 1] * temp;
                d = d * temp;
            } else {
                z[(i4 - 2 * pp) - 1] = z[(i4 + 1) - 1] * (z[(i4 - 1) - 1] / z[(i4 - 2 * pp - 2) - 1]);
                d = z[(i4 + 1) - 1] * (d / z[(i4 - 2 * pp - 2) - 1]);
            }
            emin = min(emin, &z[(i4 - 2 * pp) - 1]);
        }
        z[(4 * n0 - pp - 2) - 1] = d;
        //
        //        Now find qmax.
        //
        qmax = z[(4 * i0 - pp - 2) - 1];
        for (i4 = 4 * i0 - pp + 2; i4 <= 4 * n0 - pp - 2; i4 = i4 + 4) {
            qmax = max(qmax, &z[i4 - 1]);
        }
        //
        //        Prepare for the next iteration on K.
        //
        pp = 1 - pp;
    }
    //
    //     Initialise variables to pass to Rlasq3.
    //
    ttype = 0;
    dmin1 = zero;
    dmin2 = zero;
    dn = zero;
    dn1 = zero;
    dn2 = zero;
    g = zero;
    tau = zero;
    //
    iter = 2;
    nfail = 0;
    ndiv = 2 * (n0 - i0);
    //
    for (iwhila = 1; iwhila <= n + 1; iwhila = iwhila + 1) {
        if (n0 < 1) {
            goto statement_170;
        }
        //
        //        While array unfinished do
        //
        //        E(N0) holds the value of SIGMA when submatrix in I0:N0
        //        splits from the rest of the array, but is negated.
        //
        desig = zero;
        if (n0 == n) {
            sigma = zero;
        } else {
            sigma = -z[(4 * n0 - 1) - 1];
        }
        if (sigma < zero) {
            info = 1;
            return;
        }
        //
        //        Find last unreduced submatrix's top index I0, find QMAX and
        //        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
        //
        emax = zero;
        if (n0 > i0) {
            emin = abs(z[(4 * n0 - 5) - 1]);
        } else {
            emin = zero;
        }
        qmin = z[(4 * n0 - 3) - 1];
        qmax = qmin;
        for (i4 = 4 * n0; i4 >= 8; i4 = i4 - 4) {
            if (z[(i4 - 5) - 1] <= zero) {
                goto statement_100;
            }
            if (qmin >= four * emax) {
                qmin = min(qmin, &z[(i4 - 3) - 1]);
                emax = max(emax, &z[(i4 - 5) - 1]);
            }
            qmax = max(qmax, &z[(i4 - 7) - 1] + z[(i4 - 5) - 1]);
            emin = min(emin, &z[(i4 - 5) - 1]);
        }
        i4 = 4;
    //
    statement_100:
        i0 = i4 / 4;
        pp = 0;
        //
        if (n0 - i0 > 1) {
            dee = z[(4 * i0 - 3) - 1];
            deemin = dee;
            kmin = i0;
            for (i4 = 4 * i0 + 1; i4 <= 4 * n0 - 3; i4 = i4 + 4) {
                dee = z[i4 - 1] * (dee / (dee + z[(i4 - 2) - 1]));
                if (dee <= deemin) {
                    deemin = dee;
                    kmin = (i4 + 3) / 4;
                }
            }
            if ((kmin - i0) * 2 < n0 - kmin && deemin <= half * z[(4 * n0 - 3) - 1]) {
                ipn4 = 4 * (i0 + n0);
                pp = 2;
                for (i4 = 4 * i0; i4 <= 2 * (i0 + n0 - 1); i4 = i4 + 4) {
                    temp = z[(i4 - 3) - 1];
                    z[(i4 - 3) - 1] = z[(ipn4 - i4 - 3) - 1];
                    z[(ipn4 - i4 - 3) - 1] = temp;
                    temp = z[(i4 - 2) - 1];
                    z[(i4 - 2) - 1] = z[(ipn4 - i4 - 2) - 1];
                    z[(ipn4 - i4 - 2) - 1] = temp;
                    temp = z[(i4 - 1) - 1];
                    z[(i4 - 1) - 1] = z[(ipn4 - i4 - 5) - 1];
                    z[(ipn4 - i4 - 5) - 1] = temp;
                    temp = z[i4 - 1];
                    z[i4 - 1] = z[(ipn4 - i4 - 4) - 1];
                    z[(ipn4 - i4 - 4) - 1] = temp;
                }
            }
        }
        //
        //        Put -(initial shift) into DMIN.
        //
        dmin = -max(zero, qmin - two * sqrt(qmin) * sqrt(emax));
        //
        //        Now I0:N0 is unreduced.
        //        PP = 0 for ping, PP = 1 for pong.
        //        PP = 2 indicates that flipping was applied to the Z array and
        //               and that the tests for deflation upon entry in Rlasq3
        //               should not be performed.
        //
        nbig = 100 * (n0 - i0 + 1);
        for (iwhilb = 1; iwhilb <= nbig; iwhilb = iwhilb + 1) {
            if (i0 > n0) {
                goto statement_150;
            }
            //
            //           While submatrix unfinished take a good dqds step.
            //
            Rlasq3(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau);
            //
            pp = 1 - pp;
            //
            //           When EMIN is very small check for splits.
            //
            if (pp == 0 && n0 - i0 >= 3) {
                if (z[(4 * n0) - 1] <= tol2 * qmax || z[(4 * n0 - 1) - 1] <= tol2 * sigma) {
                    splt = i0 - 1;
                    qmax = z[(4 * i0 - 3) - 1];
                    emin = z[(4 * i0 - 1) - 1];
                    oldemn = z[(4 * i0) - 1];
                    for (i4 = 4 * i0; i4 <= 4 * (n0 - 3); i4 = i4 + 4) {
                        if (z[i4 - 1] <= tol2 * z[(i4 - 3) - 1] || z[(i4 - 1) - 1] <= tol2 * sigma) {
                            z[(i4 - 1) - 1] = -sigma;
                            splt = i4 / 4;
                            qmax = zero;
                            emin = z[(i4 + 3) - 1];
                            oldemn = z[(i4 + 4) - 1];
                        } else {
                            qmax = max(qmax, &z[(i4 + 1) - 1]);
                            emin = min(emin, &z[(i4 - 1) - 1]);
                            oldemn = min(oldemn, &z[i4 - 1]);
                        }
                    }
                    z[(4 * n0 - 1) - 1] = emin;
                    z[(4 * n0) - 1] = oldemn;
                    i0 = splt + 1;
                }
            }
            //
        }
        //
        info = 2;
        //
        //        Maximum number of iterations exceeded, restore the shift
        //        SIGMA and place the new d's and e's in a qd array.
        //        This might need to be done for several blocks
        //
        i1 = i0;
        n1 = n0;
    statement_145:
        tempq = z[(4 * i0 - 3) - 1];
        z[(4 * i0 - 3) - 1] += sigma;
        for (k = i0 + 1; k <= n0; k = k + 1) {
            tempe = z[(4 * k - 5) - 1];
            z[(4 * k - 5) - 1] = z[(4 * k - 5) - 1] * (tempq / z[(4 * k - 7) - 1]);
            tempq = z[(4 * k - 3) - 1];
            z[(4 * k - 3) - 1] += sigma + tempe - z[(4 * k - 5) - 1];
        }
        //
        //        Prepare to do this on the previous block if there is one
        //
        if (i1 > 1) {
            n1 = i1 - 1;
            while ((i1 >= 2) && (z[(4 * i1 - 5) - 1] >= zero)) {
                i1 = i1 - 1;
            }
            sigma = -z[(4 * n1 - 1) - 1];
            goto statement_145;
        }
        //
        for (k = 1; k <= n; k = k + 1) {
            z[(2 * k - 1) - 1] = z[(4 * k - 3) - 1];
            //
            //        Only the block 1..N0 is unfinished.  The rest of the e's
            //        must be essentially zero, although sometimes other data
            //        has been stored in them.
            //
            if (k < n0) {
                z[(2 * k) - 1] = z[(4 * k - 1) - 1];
            } else {
                z[(2 * k) - 1] = 0;
            }
        }
        return;
    //
    //        end IWHILB
    //
    statement_150:;
        //
    }
    //
    info = 3;
    return;
//
//     end IWHILA
//
statement_170:
    //
    //     Move q's to the front.
    //
    for (k = 2; k <= n; k = k + 1) {
        z[k - 1] = z[(4 * k - 3) - 1];
    }
    //
    //     Sort and compute sum of eigenvalues.
    //
    Rlasrt("D", n, z, iinfo);
    //
    e = zero;
    for (k = n; k >= 1; k = k - 1) {
        e += z[k - 1];
    }
    //
    //     Store trace, sum(eigenvalues) and information on performance.
    //
    z[(2 * n + 1) - 1] = trace;
    z[(2 * n + 2) - 1] = e;
    z[(2 * n + 3) - 1] = iter.real();
    z[(2 * n + 4) - 1] = ndiv.real() / pow2(n).real();
    z[(2 * n + 5) - 1] = hundrd * nfail / iter.real();
    //
    //     End of Rlasq2
    //
}
