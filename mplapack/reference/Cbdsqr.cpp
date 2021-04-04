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

void Cbdsqr(const char *uplo, INTEGER const &n, INTEGER const &ncvt, INTEGER const &nru, INTEGER const &ncc, REAL *d, REAL *e, COMPLEX *vt, INTEGER const &ldvt, COMPLEX *u, INTEGER const &ldu, COMPLEX *c, INTEGER const &ldc, REAL *rwork, INTEGER &info) {
    bool lower = false;
    bool rotate = false;
    INTEGER nm1 = 0;
    INTEGER nm12 = 0;
    INTEGER nm13 = 0;
    INTEGER idir = 0;
    REAL eps = 0.0;
    REAL unfl = 0.0;
    INTEGER i = 0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL r = 0.0;
    const REAL ten = 10.0;
    const REAL hndrd = 100.0;
    const REAL meigth = -0.125e0;
    REAL tolmul = 0.0;
    REAL tol = 0.0;
    const REAL zero = 0.0;
    REAL smax = 0.0;
    REAL sminl = 0.0;
    REAL sminoa = 0.0;
    REAL mu = 0.0;
    const INTEGER maxitr = 6;
    REAL thresh = 0.0;
    INTEGER maxit = 0;
    INTEGER iter = 0;
    INTEGER oldll = 0;
    INTEGER oldm = 0;
    INTEGER m = 0;
    REAL smin = 0.0;
    INTEGER lll = 0;
    INTEGER ll = 0;
    REAL abss = 0.0;
    REAL abse = 0.0;
    REAL sigmn = 0.0;
    REAL sigmx = 0.0;
    REAL sinr = 0.0;
    REAL cosr = 0.0;
    REAL sinl = 0.0;
    REAL cosl = 0.0;
    const REAL hndrth = 0.01e0;
    REAL shift = 0.0;
    REAL sll = 0.0;
    const REAL one = 1.0;
    REAL oldcs = 0.0;
    REAL oldsn = 0.0;
    REAL h = 0.0;
    REAL f = 0.0;
    REAL g = 0.0;
    const REAL negone = -1.0;
    INTEGER isub = 0;
    INTEGER j = 0;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    lower = Mlsame(uplo, "L");
    if (!Mlsame(uplo, "U") && !lower) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (ncvt < 0) {
        info = -3;
    } else if (nru < 0) {
        info = -4;
    } else if (ncc < 0) {
        info = -5;
    } else if ((ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < max((INTEGER)1, n))) {
        info = -9;
    } else if (ldu < max((INTEGER)1, nru)) {
        info = -11;
    } else if ((ncc == 0 && ldc < 1) || (ncc > 0 && ldc < max((INTEGER)1, n))) {
        info = -13;
    }
    if (info != 0) {
        Mxerbla("Cbdsqr", -info);
        return;
    }
    if (n == 0) {
        return;
    }
    if (n == 1) {
        goto statement_160;
    }
    //
    //     ROTATE is true if any singular vectors desired, false otherwise
    //
    rotate = (ncvt > 0) || (nru > 0) || (ncc > 0);
    //
    //     If no singular vectors desired, use qd algorithm
    //
    if (!rotate) {
        Rlasq1(n, d, e, rwork, info);
        //
        //     If INFO equals 2, dqds didn't finish, try to finish
        //
        if (info != 2) {
            return;
        }
        info = 0;
    }
    //
    nm1 = n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;
    //
    //     Get machine constants
    //
    eps = dlamch("Epsilon");
    unfl = dlamch("Safe minimum");
    //
    //     If matrix lower bidiagonal, rotate to be upper bidiagonal
    //     by applying Givens rotations on the left
    //
    if (lower) {
        for (i = 1; i <= n - 1; i = i + 1) {
            Rlartg(d[i - 1], e[i - 1], cs, sn, r);
            d[i - 1] = r;
            e[i - 1] = sn * d[(i + 1) - 1];
            d[(i + 1) - 1] = cs * d[(i + 1) - 1];
            rwork[i - 1] = cs;
            rwork[(nm1 + i) - 1] = sn;
        }
        //
        //        Update singular vectors if desired
        //
        if (nru > 0) {
            Clasr("R", "V", "F", nru, n, rwork[1 - 1], rwork[n - 1], u, ldu);
        }
        if (ncc > 0) {
            Clasr("L", "V", "F", n, ncc, rwork[1 - 1], rwork[n - 1], c, ldc);
        }
    }
    //
    //     Compute singular values to relative accuracy TOL
    //     (By setting TOL to be negative, algorithm will compute
    //     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
    //
    tolmul = max(ten, min(hndrd, pow(eps, meigth)));
    tol = tolmul * eps;
    //
    //     Compute approximate maximum, minimum singular values
    //
    smax = zero;
    for (i = 1; i <= n; i = i + 1) {
        smax = max(smax, abs(d[i - 1]));
    }
    for (i = 1; i <= n - 1; i = i + 1) {
        smax = max(smax, abs(e[i - 1]));
    }
    sminl = zero;
    if (tol >= zero) {
        //
        //        Relative accuracy desired
        //
        sminoa = abs(d[1 - 1]);
        if (sminoa == zero) {
            goto statement_50;
        }
        mu = sminoa;
        for (i = 2; i <= n; i = i + 1) {
            mu = abs(d[i - 1]) * (mu / (mu + abs(e[(i - 1) - 1])));
            sminoa = min(sminoa, mu);
            if (sminoa == zero) {
                goto statement_50;
            }
        }
    statement_50:
        sminoa = sminoa / sqrt(n.real());
        thresh = max(tol * sminoa, maxitr * n * n * unfl);
    } else {
        //
        //        Absolute accuracy desired
        //
        thresh = max(abs(tol) * smax, maxitr * n * n * unfl);
    }
    //
    //     Prepare for main iteration loop for the singular values
    //     (MAXIT is the maximum number of passes through the inner
    //     loop permitted before nonconvergence signalled.)
    //
    maxit = maxitr * n * n;
    iter = 0;
    oldll = -1;
    oldm = -1;
    //
    //     M poINTEGERs to last element of unconverged part of matrix
    //
    m = n;
//
//     Begin main iteration loop
//
statement_60:
    //
    //     Check for convergence or exceeding iteration count
    //
    if (m <= 1) {
        goto statement_160;
    }
    if (iter > maxit) {
        goto statement_200;
    }
    //
    //     Find diagonal block of matrix to work on
    //
    if (tol < zero && abs(d[m - 1]) <= thresh) {
        d[m - 1] = zero;
    }
    smax = abs(d[m - 1]);
    smin = smax;
    for (lll = 1; lll <= m - 1; lll = lll + 1) {
        ll = m - lll;
        abss = abs(d[ll - 1]);
        abse = abs(e[ll - 1]);
        if (tol < zero && abss <= thresh) {
            d[ll - 1] = zero;
        }
        if (abse <= thresh) {
            goto statement_80;
        }
        smin = min(smin, abss);
        smax = max(smax, abss, abse);
    }
    ll = 0;
    goto statement_90;
statement_80:
    e[ll - 1] = zero;
    //
    //     Matrix splits since E(LL) = 0
    //
    if (ll == m - 1) {
        //
        //        Convergence of bottom singular value, return to top of loop
        //
        m = m - 1;
        goto statement_60;
    }
statement_90:
    ll++;
    //
    //     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
    //
    if (ll == m - 1) {
        //
        //        2 by 2 block, handle separately
        //
        Rlasv2(d[(m - 1) - 1], e[(m - 1) - 1], d[m - 1], sigmn, sigmx, sinr, cosr, sinl, cosl);
        d[(m - 1) - 1] = sigmx;
        e[(m - 1) - 1] = zero;
        d[m - 1] = sigmn;
        //
        //        Compute singular vectors, if desired
        //
        if (ncvt > 0) {
            CRrot(ncvt, vt[((m - 1) - 1)], ldvt, vt[(m - 1)], ldvt, cosr, sinr);
        }
        if (nru > 0) {
            CRrot(nru, u[((m - 1) - 1) * ldu], 1, u[(m - 1) * ldu], 1, cosl, sinl);
        }
        if (ncc > 0) {
            CRrot(ncc, c[((m - 1) - 1)], ldc, c[(m - 1)], ldc, cosl, sinl);
        }
        m = m - 2;
        goto statement_60;
    }
    //
    //     If working on new submatrix, choose shift direction
    //     (from larger end diagonal element towards smaller)
    //
    if (ll > oldm || m < oldll) {
        if (abs(d[ll - 1]) >= abs(d[m - 1])) {
            //
            //           Chase bulge from top (big end) to bottom (small end)
            //
            idir = 1;
        } else {
            //
            //           Chase bulge from bottom (big end) to top (small end)
            //
            idir = 2;
        }
    }
    //
    //     Apply convergence tests
    //
    if (idir == 1) {
        //
        //        Run convergence test in forward direction
        //        First apply standard test to bottom of matrix
        //
        if (abs(e[(m - 1) - 1]) <= abs(tol) * abs(d[m - 1]) || (tol < zero && abs(e[(m - 1) - 1]) <= thresh)) {
            e[(m - 1) - 1] = zero;
            goto statement_60;
        }
        //
        if (tol >= zero) {
            //
            //           If relative accuracy desired,
            //           apply convergence criterion forward
            //
            mu = abs(d[ll - 1]);
            sminl = mu;
            for (lll = ll; lll <= m - 1; lll = lll + 1) {
                if (abs(e[lll - 1]) <= tol * mu) {
                    e[lll - 1] = zero;
                    goto statement_60;
                }
                mu = abs(d[(lll + 1) - 1]) * (mu / (mu + abs(e[lll - 1])));
                sminl = min(sminl, mu);
            }
        }
        //
    } else {
        //
        //        Run convergence test in backward direction
        //        First apply standard test to top of matrix
        //
        if (abs(e[ll - 1]) <= abs(tol) * abs(d[ll - 1]) || (tol < zero && abs(e[ll - 1]) <= thresh)) {
            e[ll - 1] = zero;
            goto statement_60;
        }
        //
        if (tol >= zero) {
            //
            //           If relative accuracy desired,
            //           apply convergence criterion backward
            //
            mu = abs(d[m - 1]);
            sminl = mu;
            for (lll = m - 1; lll >= ll; lll = lll - 1) {
                if (abs(e[lll - 1]) <= tol * mu) {
                    e[lll - 1] = zero;
                    goto statement_60;
                }
                mu = abs(d[lll - 1]) * (mu / (mu + abs(e[lll - 1])));
                sminl = min(sminl, mu);
            }
        }
    }
    oldll = ll;
    oldm = m;
    //
    //     Compute shift.  First, test if shifting would ruin relative
    //     accuracy, and if so set the shift to zero.
    //
    if (tol >= zero && n * tol * (sminl / smax) <= max(eps, hndrth * tol)) {
        //
        //        Use a zero shift to avoid loss of relative accuracy
        //
        shift = zero;
    } else {
        //
        //        Compute the shift from 2-by-2 block at end of matrix
        //
        if (idir == 1) {
            sll = abs(d[ll - 1]);
            Rlas2(d[(m - 1) - 1], e[(m - 1) - 1], d[m - 1], shift, r);
        } else {
            sll = abs(d[m - 1]);
            Rlas2(d[ll - 1], e[ll - 1], d[(ll + 1) - 1], shift, r);
        }
        //
        //        Test if shift negligible, and if so set to zero
        //
        if (sll > zero) {
            if (pow2((shift / sll)) < eps) {
                shift = zero;
            }
        }
    }
    //
    //     Increment iteration count
    //
    iter += m - ll;
    //
    //     If SHIFT = 0, do simplified QR iteration
    //
    if (shift == zero) {
        if (idir == 1) {
            //
            //           Chase bulge from top to bottom
            //           Save cosines and sines for later singular vector updates
            //
            cs = one;
            oldcs = one;
            for (i = ll; i <= m - 1; i = i + 1) {
                Rlartg(d[i - 1] * cs, e[i - 1], cs, sn, r);
                if (i > ll) {
                    e[(i - 1) - 1] = oldsn * r;
                }
                Rlartg(oldcs * r, d[(i + 1) - 1] * sn, oldcs, oldsn, d[i - 1]);
                rwork[(i - ll + 1) - 1] = cs;
                rwork[(i - ll + 1 + nm1) - 1] = sn;
                rwork[(i - ll + 1 + nm12) - 1] = oldcs;
                rwork[(i - ll + 1 + nm13) - 1] = oldsn;
            }
            h = d[m - 1] * cs;
            d[m - 1] = h * oldcs;
            e[(m - 1) - 1] = h * oldsn;
            //
            //           Update singular vectors
            //
            if (ncvt > 0) {
                Clasr("L", "V", "F", m - ll + 1, ncvt, rwork[1 - 1], rwork[n - 1], vt[(ll - 1)], ldvt);
            }
            if (nru > 0) {
                Clasr("R", "V", "F", nru, m - ll + 1, rwork[(nm12 + 1) - 1], rwork[(nm13 + 1) - 1], u[(ll - 1) * ldu], ldu);
            }
            if (ncc > 0) {
                Clasr("L", "V", "F", m - ll + 1, ncc, rwork[(nm12 + 1) - 1], rwork[(nm13 + 1) - 1], c[(ll - 1)], ldc);
            }
            //
            //           Test convergence
            //
            if (abs(e[(m - 1) - 1]) <= thresh) {
                e[(m - 1) - 1] = zero;
            }
            //
        } else {
            //
            //           Chase bulge from bottom to top
            //           Save cosines and sines for later singular vector updates
            //
            cs = one;
            oldcs = one;
            for (i = m; i >= ll + 1; i = i - 1) {
                Rlartg(d[i - 1] * cs, e[(i - 1) - 1], cs, sn, r);
                if (i < m) {
                    e[i - 1] = oldsn * r;
                }
                Rlartg(oldcs * r, d[(i - 1) - 1] * sn, oldcs, oldsn, d[i - 1]);
                rwork[(i - ll) - 1] = cs;
                rwork[(i - ll + nm1) - 1] = -sn;
                rwork[(i - ll + nm12) - 1] = oldcs;
                rwork[(i - ll + nm13) - 1] = -oldsn;
            }
            h = d[ll - 1] * cs;
            d[ll - 1] = h * oldcs;
            e[ll - 1] = h * oldsn;
            //
            //           Update singular vectors
            //
            if (ncvt > 0) {
                Clasr("L", "V", "B", m - ll + 1, ncvt, rwork[(nm12 + 1) - 1], rwork[(nm13 + 1) - 1], vt[(ll - 1)], ldvt);
            }
            if (nru > 0) {
                Clasr("R", "V", "B", nru, m - ll + 1, rwork[1 - 1], rwork[n - 1], u[(ll - 1) * ldu], ldu);
            }
            if (ncc > 0) {
                Clasr("L", "V", "B", m - ll + 1, ncc, rwork[1 - 1], rwork[n - 1], c[(ll - 1)], ldc);
            }
            //
            //           Test convergence
            //
            if (abs(e[ll - 1]) <= thresh) {
                e[ll - 1] = zero;
            }
        }
    } else {
        //
        //        Use nonzero shift
        //
        if (idir == 1) {
            //
            //           Chase bulge from top to bottom
            //           Save cosines and sines for later singular vector updates
            //
            f = (abs(d[ll - 1]) - shift) * (sign[(one - 1) + (d[ll - 1] - 1) * ldsign] + shift / d[ll - 1]);
            g = e[ll - 1];
            for (i = ll; i <= m - 1; i = i + 1) {
                Rlartg(f, g, cosr, sinr, r);
                if (i > ll) {
                    e[(i - 1) - 1] = r;
                }
                f = cosr * d[i - 1] + sinr * e[i - 1];
                e[i - 1] = cosr * e[i - 1] - sinr * d[i - 1];
                g = sinr * d[(i + 1) - 1];
                d[(i + 1) - 1] = cosr * d[(i + 1) - 1];
                Rlartg(f, g, cosl, sinl, r);
                d[i - 1] = r;
                f = cosl * e[i - 1] + sinl * d[(i + 1) - 1];
                d[(i + 1) - 1] = cosl * d[(i + 1) - 1] - sinl * e[i - 1];
                if (i < m - 1) {
                    g = sinl * e[(i + 1) - 1];
                    e[(i + 1) - 1] = cosl * e[(i + 1) - 1];
                }
                rwork[(i - ll + 1) - 1] = cosr;
                rwork[(i - ll + 1 + nm1) - 1] = sinr;
                rwork[(i - ll + 1 + nm12) - 1] = cosl;
                rwork[(i - ll + 1 + nm13) - 1] = sinl;
            }
            e[(m - 1) - 1] = f;
            //
            //           Update singular vectors
            //
            if (ncvt > 0) {
                Clasr("L", "V", "F", m - ll + 1, ncvt, rwork[1 - 1], rwork[n - 1], vt[(ll - 1)], ldvt);
            }
            if (nru > 0) {
                Clasr("R", "V", "F", nru, m - ll + 1, rwork[(nm12 + 1) - 1], rwork[(nm13 + 1) - 1], u[(ll - 1) * ldu], ldu);
            }
            if (ncc > 0) {
                Clasr("L", "V", "F", m - ll + 1, ncc, rwork[(nm12 + 1) - 1], rwork[(nm13 + 1) - 1], c[(ll - 1)], ldc);
            }
            //
            //           Test convergence
            //
            if (abs(e[(m - 1) - 1]) <= thresh) {
                e[(m - 1) - 1] = zero;
            }
            //
        } else {
            //
            //           Chase bulge from bottom to top
            //           Save cosines and sines for later singular vector updates
            //
            f = (abs(d[m - 1]) - shift) * (sign[(one - 1) + (d[m - 1] - 1) * ldsign] + shift / d[m - 1]);
            g = e[(m - 1) - 1];
            for (i = m; i >= ll + 1; i = i - 1) {
                Rlartg(f, g, cosr, sinr, r);
                if (i < m) {
                    e[i - 1] = r;
                }
                f = cosr * d[i - 1] + sinr * e[(i - 1) - 1];
                e[(i - 1) - 1] = cosr * e[(i - 1) - 1] - sinr * d[i - 1];
                g = sinr * d[(i - 1) - 1];
                d[(i - 1) - 1] = cosr * d[(i - 1) - 1];
                Rlartg(f, g, cosl, sinl, r);
                d[i - 1] = r;
                f = cosl * e[(i - 1) - 1] + sinl * d[(i - 1) - 1];
                d[(i - 1) - 1] = cosl * d[(i - 1) - 1] - sinl * e[(i - 1) - 1];
                if (i > ll + 1) {
                    g = sinl * e[(i - 2) - 1];
                    e[(i - 2) - 1] = cosl * e[(i - 2) - 1];
                }
                rwork[(i - ll) - 1] = cosr;
                rwork[(i - ll + nm1) - 1] = -sinr;
                rwork[(i - ll + nm12) - 1] = cosl;
                rwork[(i - ll + nm13) - 1] = -sinl;
            }
            e[ll - 1] = f;
            //
            //           Test convergence
            //
            if (abs(e[ll - 1]) <= thresh) {
                e[ll - 1] = zero;
            }
            //
            //           Update singular vectors if desired
            //
            if (ncvt > 0) {
                Clasr("L", "V", "B", m - ll + 1, ncvt, rwork[(nm12 + 1) - 1], rwork[(nm13 + 1) - 1], vt[(ll - 1)], ldvt);
            }
            if (nru > 0) {
                Clasr("R", "V", "B", nru, m - ll + 1, rwork[1 - 1], rwork[n - 1], u[(ll - 1) * ldu], ldu);
            }
            if (ncc > 0) {
                Clasr("L", "V", "B", m - ll + 1, ncc, rwork[1 - 1], rwork[n - 1], c[(ll - 1)], ldc);
            }
        }
    }
    //
    //     QR iteration finished, go back and check convergence
    //
    goto statement_60;
//
//     All singular values converged, so make them positive
//
statement_160:
    for (i = 1; i <= n; i = i + 1) {
        if (d[i - 1] < zero) {
            d[i - 1] = -d[i - 1];
            //
            //           Change sign of singular vectors, if desired
            //
            if (ncvt > 0) {
                CRscal(ncvt, negone, vt[(i - 1)], ldvt);
            }
        }
    }
    //
    //     Sort the singular values INTEGERo decreasing order (insertion sort on
    //     singular values, but only one transposition per singular vector)
    //
    for (i = 1; i <= n - 1; i = i + 1) {
        //
        //        Scan for smallest D(I)
        //
        isub = 1;
        smin = d[1 - 1];
        for (j = 2; j <= n + 1 - i; j = j + 1) {
            if (d[j - 1] <= smin) {
                isub = j;
                smin = d[j - 1];
            }
        }
        if (isub != n + 1 - i) {
            //
            //           Swap singular values and vectors
            //
            d[isub - 1] = d[(n + 1 - i) - 1];
            d[(n + 1 - i) - 1] = smin;
            if (ncvt > 0) {
                Cswap(ncvt, vt[(isub - 1)], ldvt, vt[((n + 1 - i) - 1)], ldvt);
            }
            if (nru > 0) {
                Cswap(nru, u[(isub - 1) * ldu], 1, u[((n + 1 - i) - 1) * ldu], 1);
            }
            if (ncc > 0) {
                Cswap(ncc, c[(isub - 1)], ldc, c[((n + 1 - i) - 1)], ldc);
            }
        }
    }
    goto statement_220;
//
//     Maximum number of iterations exceeded, failure to converge
//
statement_200:
    info = 0;
    for (i = 1; i <= n - 1; i = i + 1) {
        if (e[i - 1] != zero) {
            info++;
        }
    }
statement_220:;
    //
    //     End of Cbdsqr
    //
}
