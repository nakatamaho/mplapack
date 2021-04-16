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

void Rsbtrd(const char *vect, const char *uplo, INTEGER const n, INTEGER const kd, REAL *ab, INTEGER const ldab, REAL *d, REAL *e, REAL *q, INTEGER const ldq, REAL *work, INTEGER &info) {
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
    bool initq = Mlsame(vect, "V");
    bool wantq = initq || Mlsame(vect, "U");
    bool upper = Mlsame(uplo, "U");
    INTEGER kd1 = kd + 1;
    INTEGER kdm1 = kd - 1;
    INTEGER incx = ldab - 1;
    INTEGER iqend = 1;
    //
    info = 0;
    if (!wantq && !Mlsame(vect, "N")) {
        info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (kd < 0) {
        info = -4;
    } else if (ldab < kd1) {
        info = -6;
    } else if (ldq < max((INTEGER)1, n) && wantq) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Rsbtrd", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Initialize Q to the unit matrix, if needed
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (initq) {
        Rlaset("Full", n, n, zero, one, q, ldq);
    }
    //
    //     Wherever possible, plane rotations are generated and applied in
    //     vector operations of length NR over the index set J1:J2:KD1.
    //
    //     The cosines and sines of the plane rotations are stored in the
    //     arrays D and WORK.
    //
    INTEGER inca = kd1 * ldab;
    INTEGER kdn = min(n - 1, kd);
    INTEGER nr = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    INTEGER i = 0;
    INTEGER k = 0;
    INTEGER l = 0;
    INTEGER jend = 0;
    INTEGER jinc = 0;
    REAL temp = 0.0;
    INTEGER nrt = 0;
    INTEGER j1end = 0;
    INTEGER jin = 0;
    INTEGER lend = 0;
    INTEGER last = 0;
    INTEGER i2 = 0;
    INTEGER iqaend = 0;
    INTEGER j = 0;
    INTEGER ibl = 0;
    INTEGER iqb = 0;
    INTEGER nq = 0;
    INTEGER j1inc = 0;
    if (upper) {
        //
        if (kd > 1) {
            //
            //           Reduce to tridiagonal form, working with upper triangle
            //
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;
            //
            for (i = 1; i <= n - 2; i = i + 1) {
                //
                //              Reduce i-th row of matrix to tridiagonal form
                //
                for (k = kdn + 1; k >= 2; k = k - 1) {
                    j1 += kdn;
                    j2 += kdn;
                    //
                    if (nr > 0) {
                        //
                        //                    generate plane rotations to annihilate nonzero
                        //                    elements which have been created outside the band
                        //
                        Rlargv(nr, &ab[((j1 - 1) - 1) * ldab], inca, &work[j1 - 1], kd1, &d[j1 - 1], kd1);
                        //
                        //                    apply rotations from the right
                        //
                        //                    Dependent on the the number of diagonals either
                        //                    Rlartv or Rrot is used
                        //
                        if (nr >= 2 * kd - 1) {
                            for (l = 1; l <= kd - 1; l = l + 1) {
                                Rlartv(nr, &ab[((l + 1) - 1) + ((j1 - 1) - 1) * ldab], inca, &ab[(l - 1) + (j1 - 1) * ldab], inca, &d[j1 - 1], &work[j1 - 1], kd1);
                            }
                            //
                        } else {
                            jend = j1 + (nr - 1) * kd1;
                            for (jinc = j1; jinc <= jend; jinc = jinc + kd1) {
                                Rrot(kdm1, &ab[(2 - 1) + ((jinc - 1) - 1) * ldab], 1, &ab[(jinc - 1) * ldab], 1, d[jinc - 1], work[jinc - 1]);
                            }
                        }
                    }
                    //
                    if (k > 2) {
                        if (k <= n - i + 1) {
                            //
                            //                       generate plane rotation to annihilate a(i,i+k-1)
                            //                       within the band
                            //
                            Rlartg(ab[((kd - k + 3) - 1) + ((i + k - 2) - 1) * ldab], ab[((kd - k + 2) - 1) + ((i + k - 1) - 1) * ldab], d[(i + k - 1) - 1], work[(i + k - 1) - 1], temp);
                            ab[((kd - k + 3) - 1) + ((i + k - 2) - 1) * ldab] = temp;
                            //
                            //                       apply rotation from the right
                            //
                            Rrot(k - 3, &ab[((kd - k + 4) - 1) + ((i + k - 2) - 1) * ldab], 1, &ab[((kd - k + 3) - 1) + ((i + k - 1) - 1) * ldab], 1, d[(i + k - 1) - 1], work[(i + k - 1) - 1]);
                        }
                        nr++;
                        j1 = j1 - kdn - 1;
                    }
                    //
                    //                 apply plane rotations from both sides to diagonal
                    //                 blocks
                    //
                    if (nr > 0) {
                        Rlar2v(nr, &ab[(kd1 - 1) + ((j1 - 1) - 1) * ldab], &ab[(kd1 - 1) + (j1 - 1) * ldab], &ab[(kd - 1) + (j1 - 1) * ldab], inca, &d[j1 - 1], &work[j1 - 1], kd1);
                    }
                    //
                    //                 apply plane rotations from the left
                    //
                    if (nr > 0) {
                        if (2 * kd - 1 < nr) {
                            //
                            //                    Dependent on the the number of diagonals either
                            //                    Rlartv or Rrot is used
                            //
                            for (l = 1; l <= kd - 1; l = l + 1) {
                                if (j2 + l > n) {
                                    nrt = nr - 1;
                                } else {
                                    nrt = nr;
                                }
                                if (nrt > 0) {
                                    Rlartv(nrt, &ab[((kd - l) - 1) + ((j1 + l) - 1) * ldab], inca, &ab[((kd - l + 1) - 1) + ((j1 + l) - 1) * ldab], inca, &d[j1 - 1], &work[j1 - 1], kd1);
                                }
                            }
                        } else {
                            j1end = j1 + kd1 * (nr - 2);
                            if (j1end >= j1) {
                                for (jin = j1; jin <= j1end; jin = jin + kd1) {
                                    Rrot(kd - 1, &ab[((kd - 1) - 1) + ((jin + 1) - 1) * ldab], incx, &ab[(kd - 1) + ((jin + 1) - 1) * ldab], incx, d[jin - 1], work[jin - 1]);
                                }
                            }
                            lend = min(kdm1, n - j2);
                            last = j1end + kd1;
                            if (lend > 0) {
                                Rrot(lend, &ab[((kd - 1) - 1) + ((last + 1) - 1) * ldab], incx, &ab[(kd - 1) + ((last + 1) - 1) * ldab], incx, d[last - 1], work[last - 1]);
                            }
                        }
                    }
                    //
                    if (wantq) {
                        //
                        //                    accumulate product of plane rotations in Q
                        //
                        if (initq) {
                            //
                            //                 take advantage of the fact that Q was
                            //                 initially the Identity matrix
                            //
                            iqend = max(iqend, j2);
                            i2 = max(0, k - 3);
                            iqaend = 1 + i * kd;
                            if (k == 2) {
                                iqaend += kd;
                            }
                            iqaend = min(iqaend, iqend);
                            for (j = j1; j <= j2; j = j + kd1) {
                                ibl = i - i2 / kdm1;
                                i2++;
                                iqb = max((INTEGER)1, j - ibl);
                                nq = 1 + iqaend - iqb;
                                iqaend = min(iqaend + kd, iqend);
                                Rrot(nq, &q[(iqb - 1) + ((j - 1) - 1) * ldq], 1, &q[(iqb - 1) + (j - 1) * ldq], 1, d[j - 1], work[j - 1]);
                            }
                        } else {
                            //
                            for (j = j1; j <= j2; j = j + kd1) {
                                Rrot(n, &q[((j - 1) - 1) * ldq], 1, &q[(j - 1) * ldq], 1, d[j - 1], work[j - 1]);
                            }
                        }
                        //
                    }
                    //
                    if (j2 + kdn > n) {
                        //
                        //                    adjust J2 to keep within the bounds of the matrix
                        //
                        nr = nr - 1;
                        j2 = j2 - kdn - 1;
                    }
                    //
                    for (j = j1; j <= j2; j = j + kd1) {
                        //
                        //                    create nonzero element a(j-1,j+kd) outside the band
                        //                    and store it in WORK
                        //
                        work[(j + kd) - 1] = work[j - 1] * ab[((j + kd) - 1) * ldab];
                        ab[((j + kd) - 1) * ldab] = d[j - 1] * ab[((j + kd) - 1) * ldab];
                    }
                }
            }
        }
        //
        if (kd > 0) {
            //
            //           copy off-diagonal elements to E
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                e[i - 1] = ab[(kd - 1) + ((i + 1) - 1) * ldab];
            }
        } else {
            //
            //           set E to zero if original matrix was diagonal
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                e[i - 1] = zero;
            }
        }
        //
        //        copy diagonal elements to D
        //
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = ab[(kd1 - 1) + (i - 1) * ldab];
        }
        //
    } else {
        //
        if (kd > 1) {
            //
            //           Reduce to tridiagonal form, working with lower triangle
            //
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;
            //
            for (i = 1; i <= n - 2; i = i + 1) {
                //
                //              Reduce i-th column of matrix to tridiagonal form
                //
                for (k = kdn + 1; k >= 2; k = k - 1) {
                    j1 += kdn;
                    j2 += kdn;
                    //
                    if (nr > 0) {
                        //
                        //                    generate plane rotations to annihilate nonzero
                        //                    elements which have been created outside the band
                        //
                        Rlargv(nr, &ab[(kd1 - 1) + ((j1 - kd1) - 1) * ldab], inca, &work[j1 - 1], kd1, &d[j1 - 1], kd1);
                        //
                        //                    apply plane rotations from one side
                        //
                        //                    Dependent on the the number of diagonals either
                        //                    Rlartv or Rrot is used
                        //
                        if (nr > 2 * kd - 1) {
                            for (l = 1; l <= kd - 1; l = l + 1) {
                                Rlartv(nr, &ab[((kd1 - l) - 1) + ((j1 - kd1 + l) - 1) * ldab], inca, &ab[((kd1 - l + 1) - 1) + ((j1 - kd1 + l) - 1) * ldab], inca, &d[j1 - 1], &work[j1 - 1], kd1);
                            }
                        } else {
                            jend = j1 + kd1 * (nr - 1);
                            for (jinc = j1; jinc <= jend; jinc = jinc + kd1) {
                                Rrot(kdm1, &ab[(kd - 1) + ((jinc - kd) - 1) * ldab], incx, &ab[(kd1 - 1) + ((jinc - kd) - 1) * ldab], incx, d[jinc - 1], work[jinc - 1]);
                            }
                        }
                        //
                    }
                    //
                    if (k > 2) {
                        if (k <= n - i + 1) {
                            //
                            //                       generate plane rotation to annihilate a(i+k-1,i)
                            //                       within the band
                            //
                            Rlartg(ab[((k - 1) - 1) + (i - 1) * ldab], ab[(k - 1) + (i - 1) * ldab], d[(i + k - 1) - 1], work[(i + k - 1) - 1], temp);
                            ab[((k - 1) - 1) + (i - 1) * ldab] = temp;
                            //
                            //                       apply rotation from the left
                            //
                            Rrot(k - 3, &ab[((k - 2) - 1) + ((i + 1) - 1) * ldab], ldab - 1, &ab[((k - 1) - 1) + ((i + 1) - 1) * ldab], ldab - 1, d[(i + k - 1) - 1], work[(i + k - 1) - 1]);
                        }
                        nr++;
                        j1 = j1 - kdn - 1;
                    }
                    //
                    //                 apply plane rotations from both sides to diagonal
                    //                 blocks
                    //
                    if (nr > 0) {
                        Rlar2v(nr, &ab[((j1 - 1) - 1) * ldab], &ab[(j1 - 1) * ldab], &ab[(2 - 1) + ((j1 - 1) - 1) * ldab], inca, &d[j1 - 1], &work[j1 - 1], kd1);
                    }
                    //
                    //                 apply plane rotations from the right
                    //
                    //                    Dependent on the the number of diagonals either
                    //                    Rlartv or Rrot is used
                    //
                    if (nr > 0) {
                        if (nr > 2 * kd - 1) {
                            for (l = 1; l <= kd - 1; l = l + 1) {
                                if (j2 + l > n) {
                                    nrt = nr - 1;
                                } else {
                                    nrt = nr;
                                }
                                if (nrt > 0) {
                                    Rlartv(nrt, &ab[((l + 2) - 1) + ((j1 - 1) - 1) * ldab], inca, &ab[((l + 1) - 1) + (j1 - 1) * ldab], inca, &d[j1 - 1], &work[j1 - 1], kd1);
                                }
                            }
                        } else {
                            j1end = j1 + kd1 * (nr - 2);
                            if (j1end >= j1) {
                                for (j1inc = j1; j1inc <= j1end; j1inc = j1inc + kd1) {
                                    Rrot(kdm1, &ab[(3 - 1) + ((j1inc - 1) - 1) * ldab], 1, &ab[(2 - 1) + (j1inc - 1) * ldab], 1, d[j1inc - 1], work[j1inc - 1]);
                                }
                            }
                            lend = min(kdm1, n - j2);
                            last = j1end + kd1;
                            if (lend > 0) {
                                Rrot(lend, &ab[(3 - 1) + ((last - 1) - 1) * ldab], 1, &ab[(2 - 1) + (last - 1) * ldab], 1, d[last - 1], work[last - 1]);
                            }
                        }
                    }
                    //
                    if (wantq) {
                        //
                        //                    accumulate product of plane rotations in Q
                        //
                        if (initq) {
                            //
                            //                 take advantage of the fact that Q was
                            //                 initially the Identity matrix
                            //
                            iqend = max(iqend, j2);
                            i2 = max(0, k - 3);
                            iqaend = 1 + i * kd;
                            if (k == 2) {
                                iqaend += kd;
                            }
                            iqaend = min(iqaend, iqend);
                            for (j = j1; j <= j2; j = j + kd1) {
                                ibl = i - i2 / kdm1;
                                i2++;
                                iqb = max((INTEGER)1, j - ibl);
                                nq = 1 + iqaend - iqb;
                                iqaend = min(iqaend + kd, iqend);
                                Rrot(nq, &q[(iqb - 1) + ((j - 1) - 1) * ldq], 1, &q[(iqb - 1) + (j - 1) * ldq], 1, d[j - 1], work[j - 1]);
                            }
                        } else {
                            //
                            for (j = j1; j <= j2; j = j + kd1) {
                                Rrot(n, &q[((j - 1) - 1) * ldq], 1, &q[(j - 1) * ldq], 1, d[j - 1], work[j - 1]);
                            }
                        }
                    }
                    //
                    if (j2 + kdn > n) {
                        //
                        //                    adjust J2 to keep within the bounds of the matrix
                        //
                        nr = nr - 1;
                        j2 = j2 - kdn - 1;
                    }
                    //
                    for (j = j1; j <= j2; j = j + kd1) {
                        //
                        //                    create nonzero element a(j+kd,j-1) outside the
                        //                    band and store it in WORK
                        //
                        work[(j + kd) - 1] = work[j - 1] * ab[(kd1 - 1) + (j - 1) * ldab];
                        ab[(kd1 - 1) + (j - 1) * ldab] = d[j - 1] * ab[(kd1 - 1) + (j - 1) * ldab];
                    }
                }
            }
        }
        //
        if (kd > 0) {
            //
            //           copy off-diagonal elements to E
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                e[i - 1] = ab[(2 - 1) + (i - 1) * ldab];
            }
        } else {
            //
            //           set E to zero if original matrix was diagonal
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                e[i - 1] = zero;
            }
        }
        //
        //        copy diagonal elements to D
        //
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = ab[(i - 1) * ldab];
        }
    }
    //
    //     End of Rsbtrd
    //
}
