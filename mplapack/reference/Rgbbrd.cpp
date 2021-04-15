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

void Rgbbrd(const char *vect, INTEGER const m, INTEGER const n, INTEGER const ncc, INTEGER const kl, INTEGER const ku, REAL *ab, INTEGER const ldab, REAL *d, REAL *e, REAL *q, INTEGER const ldq, REAL *pt, INTEGER const ldpt, REAL *c, INTEGER const ldc, REAL *work, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    bool wantb = Mlsame(vect, "B");
    bool wantq = Mlsame(vect, "Q") || wantb;
    bool wantpt = Mlsame(vect, "P") || wantb;
    bool wantc = ncc > 0;
    INTEGER klu1 = kl + ku + 1;
    info = 0;
    if (!wantq && !wantpt && !Mlsame(vect, "N")) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ncc < 0) {
        info = -4;
    } else if (kl < 0) {
        info = -5;
    } else if (ku < 0) {
        info = -6;
    } else if (ldab < klu1) {
        info = -8;
    } else if (ldq < 1 || wantq && ldq < max((INTEGER)1, m)) {
        info = -12;
    } else if (ldpt < 1 || wantpt && ldpt < max((INTEGER)1, n)) {
        info = -14;
    } else if (ldc < 1 || wantc && ldc < max((INTEGER)1, m)) {
        info = -16;
    }
    if (info != 0) {
        Mxerbla("Rgbbrd", -info);
        return;
    }
    //
    //     Initialize Q and P**T to the unit matrix, if needed
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (wantq) {
        Rlaset("Full", m, m, zero, one, q, ldq);
    }
    if (wantpt) {
        Rlaset("Full", n, n, zero, one, pt, ldpt);
    }
    //
    //     Quick return if possible.
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    INTEGER minmn = min(m, n);
    //
    INTEGER ml0 = 0;
    INTEGER mu0 = 0;
    INTEGER mn = 0;
    INTEGER klm = 0;
    INTEGER kun = 0;
    INTEGER kb = 0;
    INTEGER kb1 = 0;
    INTEGER inca = 0;
    INTEGER nr = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    INTEGER i = 0;
    INTEGER ml = 0;
    INTEGER mu = 0;
    INTEGER kk = 0;
    INTEGER l = 0;
    INTEGER nrt = 0;
    REAL ra = 0.0;
    INTEGER j = 0;
    if (kl + ku > 1) {
        //
        //        Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce
        //        first to lower bidiagonal form and then transform to upper
        //        bidiagonal
        //
        if (ku > 0) {
            ml0 = 1;
            mu0 = 2;
        } else {
            ml0 = 2;
            mu0 = 1;
        }
        //
        //        Wherever possible, plane rotations are generated and applied in
        //        vector operations of length NR over the index set J1:J2:KLU1.
        //
        //        The sines of the plane rotations are stored in WORK(1:max(m,n))
        //        and the cosines in WORK(max(m,n)+1:2*max(m,n)).
        //
        mn = max(m, n);
        klm = min(m - 1, kl);
        kun = min(n - 1, ku);
        kb = klm + kun;
        kb1 = kb + 1;
        inca = kb1 * ldab;
        nr = 0;
        j1 = klm + 2;
        j2 = 1 - kun;
        //
        for (i = 1; i <= minmn; i = i + 1) {
            //
            //           Reduce i-th column and i-th row of matrix to bidiagonal form
            //
            ml = klm + 1;
            mu = kun + 1;
            for (kk = 1; kk <= kb; kk = kk + 1) {
                j1 += kb;
                j2 += kb;
                //
                //              generate plane rotations to annihilate nonzero elements
                //              which have been created below the band
                //
                if (nr > 0) {
                    Rlargv(nr, ab[(klu1 - 1) + ((j1 - klm - 1) - 1) * ldab], inca, &work[j1 - 1], kb1, &work[(mn + j1) - 1], kb1);
                }
                //
                //              apply plane rotations from the left
                //
                for (l = 1; l <= kb; l = l + 1) {
                    if (j2 - klm + l - 1 > n) {
                        nrt = nr - 1;
                    } else {
                        nrt = nr;
                    }
                    if (nrt > 0) {
                        Rlartv(nrt, ab[((klu1 - l) - 1) + ((j1 - klm + l - 1) - 1) * ldab], inca, ab[((klu1 - l + 1) - 1) + ((j1 - klm + l - 1) - 1) * ldab], inca, &work[(mn + j1) - 1], &work[j1 - 1], kb1);
                    }
                }
                //
                if (ml > ml0) {
                    if (ml <= m - i + 1) {
                        //
                        //                    generate plane rotation to annihilate a(i+ml-1,i)
                        //                    within the band, and apply rotation from the left
                        //
                        Rlartg(ab[((ku + ml - 1) - 1) + (i - 1) * ldab], ab[((ku + ml) - 1) + (i - 1) * ldab], &work[(mn + i + ml - 1) - 1], &work[(i + ml - 1) - 1], ra);
                        ab[((ku + ml - 1) - 1) + (i - 1) * ldab] = ra;
                        if (i < n) {
                            Rrot(min(ku + ml - 2, n - i), ab[((ku + ml - 2) - 1) + ((i + 1) - 1) * ldab], ldab - 1, ab[((ku + ml - 1) - 1) + ((i + 1) - 1) * ldab], ldab - 1, &work[(mn + i + ml - 1) - 1], &work[(i + ml - 1) - 1]);
                        }
                    }
                    nr++;
                    j1 = j1 - kb1;
                }
                //
                if (wantq) {
                    //
                    //                 accumulate product of plane rotations in Q
                    //
                    for (j = j1; j <= j2; j = j + kb1) {
                        Rrot(m, q[((j - 1) - 1) * ldq], 1, q[(j - 1) * ldq], 1, &work[(mn + j) - 1], &work[j - 1]);
                    }
                }
                //
                if (wantc) {
                    //
                    //                 apply plane rotations to C
                    //
                    for (j = j1; j <= j2; j = j + kb1) {
                        Rrot(ncc, &c[((j - 1) - 1)], ldc, &c[(j - 1)], ldc, &work[(mn + j) - 1], &work[j - 1]);
                    }
                }
                //
                if (j2 + kun > n) {
                    //
                    //                 adjust J2 to keep within the bounds of the matrix
                    //
                    nr = nr - 1;
                    j2 = j2 - kb1;
                }
                //
                for (j = j1; j <= j2; j = j + kb1) {
                    //
                    //                 create nonzero element a(j-1,j+ku) above the band
                    //                 and store it in WORK(n+1:2*n)
                    //
                    work[(j + kun) - 1] = work[j - 1] * ab[((j + kun) - 1) * ldab];
                    ab[((j + kun) - 1) * ldab] = work[(mn + j) - 1] * ab[((j + kun) - 1) * ldab];
                }
                //
                //              generate plane rotations to annihilate nonzero elements
                //              which have been generated above the band
                //
                if (nr > 0) {
                    Rlargv(nr, ab[((j1 + kun - 1) - 1) * ldab], inca, &work[(j1 + kun) - 1], kb1, &work[(mn + j1 + kun) - 1], kb1);
                }
                //
                //              apply plane rotations from the right
                //
                for (l = 1; l <= kb; l = l + 1) {
                    if (j2 + l - 1 > m) {
                        nrt = nr - 1;
                    } else {
                        nrt = nr;
                    }
                    if (nrt > 0) {
                        Rlartv(nrt, ab[((l + 1) - 1) + ((j1 + kun - 1) - 1) * ldab], inca, ab[(l - 1) + ((j1 + kun) - 1) * ldab], inca, &work[(mn + j1 + kun) - 1], &work[(j1 + kun) - 1], kb1);
                    }
                }
                //
                if (ml == ml0 && mu > mu0) {
                    if (mu <= n - i + 1) {
                        //
                        //                    generate plane rotation to annihilate a(i,i+mu-1)
                        //                    within the band, and apply rotation from the right
                        //
                        Rlartg(ab[((ku - mu + 3) - 1) + ((i + mu - 2) - 1) * ldab], ab[((ku - mu + 2) - 1) + ((i + mu - 1) - 1) * ldab], &work[(mn + i + mu - 1) - 1], &work[(i + mu - 1) - 1], ra);
                        ab[((ku - mu + 3) - 1) + ((i + mu - 2) - 1) * ldab] = ra;
                        Rrot(min(kl + mu - 2, m - i), ab[((ku - mu + 4) - 1) + ((i + mu - 2) - 1) * ldab], 1, ab[((ku - mu + 3) - 1) + ((i + mu - 1) - 1) * ldab], 1, &work[(mn + i + mu - 1) - 1], &work[(i + mu - 1) - 1]);
                    }
                    nr++;
                    j1 = j1 - kb1;
                }
                //
                if (wantpt) {
                    //
                    //                 accumulate product of plane rotations in P**T
                    //
                    for (j = j1; j <= j2; j = j + kb1) {
                        Rrot(n, pt[((j + kun - 1) - 1)], ldpt, pt[((j + kun) - 1)], ldpt, &work[(mn + j + kun) - 1], &work[(j + kun) - 1]);
                    }
                }
                //
                if (j2 + kb > m) {
                    //
                    //                 adjust J2 to keep within the bounds of the matrix
                    //
                    nr = nr - 1;
                    j2 = j2 - kb1;
                }
                //
                for (j = j1; j <= j2; j = j + kb1) {
                    //
                    //                 create nonzero element a(j+kl+ku,j+ku-1) below the
                    //                 band and store it in WORK(1:n)
                    //
                    work[(j + kb) - 1] = work[(j + kun) - 1] * ab[(klu1 - 1) + ((j + kun) - 1) * ldab];
                    ab[(klu1 - 1) + ((j + kun) - 1) * ldab] = work[(mn + j + kun) - 1] * ab[(klu1 - 1) + ((j + kun) - 1) * ldab];
                }
                //
                if (ml > ml0) {
                    ml = ml - 1;
                } else {
                    mu = mu - 1;
                }
            }
        }
    }
    //
    REAL rc = 0.0;
    REAL rs = 0.0;
    REAL rb = 0.0;
    if (ku == 0 && kl > 0) {
        //
        //        A has been reduced to lower bidiagonal form
        //
        //        Transform lower bidiagonal form to upper bidiagonal by applying
        //        plane rotations from the left, storing diagonal elements in D
        //        and off-diagonal elements in E
        //
        for (i = 1; i <= min(m - 1, n); i = i + 1) {
            Rlartg(ab[(i - 1) * ldab], ab[(2 - 1) + (i - 1) * ldab], rc, rs, ra);
            d[i - 1] = ra;
            if (i < n) {
                e[i - 1] = rs * ab[((i + 1) - 1) * ldab];
                ab[((i + 1) - 1) * ldab] = rc * ab[((i + 1) - 1) * ldab];
            }
            if (wantq) {
                Rrot(m, q[(i - 1) * ldq], 1, q[((i + 1) - 1) * ldq], 1, rc, rs);
            }
            if (wantc) {
                Rrot(ncc, &c[(i - 1)], ldc, &c[((i + 1) - 1)], ldc, rc, rs);
            }
        }
        if (m <= n) {
            d[m - 1] = ab[(m - 1) * ldab];
        }
    } else if (ku > 0) {
        //
        //        A has been reduced to upper bidiagonal form
        //
        if (m < n) {
            //
            //           Annihilate a(m,m+1) by applying plane rotations from the
            //           right, storing diagonal elements in D and off-diagonal
            //           elements in E
            //
            rb = ab[(ku - 1) + ((m + 1) - 1) * ldab];
            for (i = m; i >= 1; i = i - 1) {
                Rlartg(ab[((ku + 1) - 1) + (i - 1) * ldab], rb, rc, rs, ra);
                d[i - 1] = ra;
                if (i > 1) {
                    rb = -rs * ab[(ku - 1) + (i - 1) * ldab];
                    e[(i - 1) - 1] = rc * ab[(ku - 1) + (i - 1) * ldab];
                }
                if (wantpt) {
                    Rrot(n, pt[(i - 1)], ldpt, pt[((m + 1) - 1)], ldpt, rc, rs);
                }
            }
        } else {
            //
            //           Copy off-diagonal elements to E and diagonal elements to D
            //
            for (i = 1; i <= minmn - 1; i = i + 1) {
                e[i - 1] = ab[(ku - 1) + ((i + 1) - 1) * ldab];
            }
            for (i = 1; i <= minmn; i = i + 1) {
                d[i - 1] = ab[((ku + 1) - 1) + (i - 1) * ldab];
            }
        }
    } else {
        //
        //        A is diagonal. Set elements of E to zero and copy diagonal
        //        elements to D.
        //
        for (i = 1; i <= minmn - 1; i = i + 1) {
            e[i - 1] = zero;
        }
        for (i = 1; i <= minmn; i = i + 1) {
            d[i - 1] = ab[(i - 1) * ldab];
        }
    }
    //
    //     End of Rgbbrd
    //
}
