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

void Cgghd3(const char *compq, const char *compz, INTEGER const &n, INTEGER const &ilo, INTEGER const &ihi, COMPLEX *a, INTEGER const &lda, COMPLEX *b, INTEGER const &ldb, COMPLEX *q, INTEGER const &ldq, COMPLEX *z, INTEGER const &ldz, COMPLEX *work, INTEGER const &lwork, INTEGER &info) {
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
    //     Decode and test the input parameters.
    //
    info = 0;
    INTEGER nb = iMlaenv[("Cgghd3" - 1) * ldiMlaenv];
    INTEGER lwkopt = max(6 * n * nb, 1);
    work[1 - 1] = COMPLEX(lwkopt);
    bool initq = Mlsame(compq, "I");
    bool wantq = initq || Mlsame(compq, "V");
    bool initz = Mlsame(compz, "I");
    bool wantz = initz || Mlsame(compz, "V");
    bool lquery = (lwork == -1);
    //
    if (!Mlsame(compq, "N") && !wantq) {
        info = -1;
    } else if (!Mlsame(compz, "N") && !wantz) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ilo < 1) {
        info = -4;
    } else if (ihi > n || ihi < ilo - 1) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if ((wantq && ldq < n) || ldq < 1) {
        info = -11;
    } else if ((wantz && ldz < n) || ldz < 1) {
        info = -13;
    } else if (lwork < 1 && !lquery) {
        info = -15;
    }
    if (info != 0) {
        Mxerbla("Cgghd3", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Initialize Q and Z if desired.
    //
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    if (initq) {
        Claset("All", n, n, czero, cone, q, ldq);
    }
    if (initz) {
        Claset("All", n, n, czero, cone, z, ldz);
    }
    //
    //     Zero out lower triangle of B.
    //
    if (n > 1) {
        Claset("Lower", n - 1, n - 1, czero, czero, b[(2 - 1)], ldb);
    }
    //
    //     Quick return if possible
    //
    INTEGER nh = ihi - ilo + 1;
    if (nh <= 1) {
        work[1 - 1] = cone;
        return;
    }
    //
    //     Determine the blocksize.
    //
    INTEGER nbmin = iMlaenv[(2 - 1) + ("Cgghd3" - 1) * ldiMlaenv];
    INTEGER nx = 0;
    if (nb > 1 && nb < nh) {
        //
        //        Determine when to use unblocked instead of blocked code.
        //
        nx = max(nb, iMlaenv[(3 - 1) + ("Cgghd3" - 1) * ldiMlaenv]);
        if (nx < nh) {
            //
            //           Determine if workspace is large enough for blocked code.
            //
            if (lwork < lwkopt) {
                //
                //              Not enough workspace to use optimal NB:  determine the
                //              minimum value of NB, and reduce NB or force use of
                //              unblocked code.
                //
                nbmin = max(2, iMlaenv[(2 - 1) + ("Cgghd3" - 1) * ldiMlaenv]);
                if (lwork >= 6 * n * nbmin) {
                    nb = lwork / (6 * n);
                } else {
                    nb = 1;
                }
            }
        }
    }
    //
    INTEGER jcol = 0;
    INTEGER kacc22 = 0;
    bool blk22 = false;
    INTEGER nnb = 0;
    INTEGER n2nb = 0;
    INTEGER nblst = 0;
    INTEGER pw = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    COMPLEX temp = 0.0;
    REAL c = 0.0;
    COMPLEX s = 0.0;
    INTEGER ppw = 0;
    INTEGER len = 0;
    INTEGER jrow = 0;
    COMPLEX ctemp = 0.0;
    INTEGER jj = 0;
    INTEGER ppwo = 0;
    INTEGER j0 = 0;
    INTEGER top = 0;
    COMPLEX c1 = 0.0;
    COMPLEX s1 = 0.0;
    COMPLEX c2 = 0.0;
    COMPLEX s2 = 0.0;
    INTEGER k = 0;
    COMPLEX temp1 = 0.0;
    COMPLEX temp2 = 0.0;
    COMPLEX temp3 = 0.0;
    INTEGER cola = 0;
    INTEGER ierr = 0;
    INTEGER topq = 0;
    if (nb < nbmin || nb >= nh) {
        //
        //        Use unblocked code below
        //
        jcol = ilo;
        //
    } else {
        //
        //        Use blocked code
        //
        kacc22 = iMlaenv[(16 - 1) + ("Cgghd3" - 1) * ldiMlaenv];
        blk22 = kacc22 == 2;
        for (jcol = ilo; jcol <= ihi - 2; jcol = jcol + nb) {
            nnb = min(nb, ihi - jcol - 1);
            //
            //           Initialize small unitary factors that will hold the
            //           accumulated Givens rotations in workspace.
            //           N2NB   denotes the number of 2*NNB-by-2*NNB factors
            //           NBLST  denotes the (possibly smaller) order of the last
            //                  factor.
            //
            n2nb = (ihi - jcol - 1) / nnb - 1;
            nblst = ihi - jcol - n2nb * nnb;
            Claset("All", nblst, nblst, czero, cone, work, nblst);
            pw = nblst * nblst + 1;
            for (i = 1; i <= n2nb; i = i + 1) {
                Claset("All", 2 * nnb, 2 * nnb, czero, cone, work[pw - 1], 2 * nnb);
                pw += 4 * nnb * nnb;
            }
            //
            //           Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form.
            //
            for (j = jcol; j <= jcol + nnb - 1; j = j + 1) {
                //
                //              Reduce Jth column of A. Store cosines and sines in Jth
                //              column of A and B, respectively.
                //
                for (i = ihi; i >= j + 2; i = i - 1) {
                    temp = a[((i - 1) - 1) + (j - 1) * lda];
                    Clartg(temp, a[(i - 1) + (j - 1) * lda], c, s, a[((i - 1) - 1) + (j - 1) * lda]);
                    a[(i - 1) + (j - 1) * lda] = COMPLEX(c);
                    b[(i - 1) + (j - 1) * ldb] = s;
                }
                //
                //              Accumulate Givens rotations INTEGERo workspace array.
                //
                ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
                len = 2 + j - jcol;
                jrow = j + n2nb * nnb + 2;
                for (i = ihi; i >= jrow; i = i - 1) {
                    ctemp = a[(i - 1) + (j - 1) * lda];
                    s = b[(i - 1) + (j - 1) * ldb];
                    for (jj = ppw; jj <= ppw + len - 1; jj = jj + 1) {
                        temp = work[(jj + nblst) - 1];
                        work[(jj + nblst) - 1] = ctemp * temp - s * work[jj - 1];
                        work[jj - 1] = conj(s) * temp + ctemp * work[jj - 1];
                    }
                    len++;
                    ppw = ppw - nblst - 1;
                }
                //
                ppwo = nblst * nblst + (nnb + j - jcol - 1) * 2 * nnb + nnb;
                j0 = jrow - nnb;
                for (jrow = j0; jrow <= j + 2; jrow = jrow + -nnb) {
                    ppw = ppwo;
                    len = 2 + j - jcol;
                    for (i = jrow + nnb - 1; i >= jrow; i = i - 1) {
                        ctemp = a[(i - 1) + (j - 1) * lda];
                        s = b[(i - 1) + (j - 1) * ldb];
                        for (jj = ppw; jj <= ppw + len - 1; jj = jj + 1) {
                            temp = work[(jj + 2 * nnb) - 1];
                            work[(jj + 2 * nnb) - 1] = ctemp * temp - s * work[jj - 1];
                            work[jj - 1] = conj(s) * temp + ctemp * work[jj - 1];
                        }
                        len++;
                        ppw = ppw - 2 * nnb - 1;
                    }
                    ppwo += 4 * nnb * nnb;
                }
                //
                //              TOP denotes the number of top rows in A and B that will
                //              not be updated during the next steps.
                //
                if (jcol <= 2) {
                    top = 0;
                } else {
                    top = jcol;
                }
                //
                //              Propagate transformations through B and replace stored
                //              left sines/cosines by right sines/cosines.
                //
                for (jj = n; jj >= j + 1; jj = jj - 1) {
                    //
                    //                 Update JJth column of B.
                    //
                    for (i = min(jj + 1, ihi); i >= j + 2; i = i - 1) {
                        ctemp = a[(i - 1) + (j - 1) * lda];
                        s = b[(i - 1) + (j - 1) * ldb];
                        temp = b[(i - 1) + (jj - 1) * ldb];
                        b[(i - 1) + (jj - 1) * ldb] = ctemp * temp - conj(s) * b[((i - 1) - 1) + (jj - 1) * ldb];
                        b[((i - 1) - 1) + (jj - 1) * ldb] = s * temp + ctemp * b[((i - 1) - 1) + (jj - 1) * ldb];
                    }
                    //
                    //                 Annihilate B( JJ+1, JJ ).
                    //
                    if (jj < ihi) {
                        temp = b[((jj + 1) - 1) + ((jj + 1) - 1) * ldb];
                        Clartg(temp, b[((jj + 1) - 1) + (jj - 1) * ldb], c, s, b[((jj + 1) - 1) + ((jj + 1) - 1) * ldb]);
                        b[((jj + 1) - 1) + (jj - 1) * ldb] = czero;
                        Crot(jj - top, b[((top + 1) - 1) + ((jj + 1) - 1) * ldb], 1, b[((top + 1) - 1) + (jj - 1) * ldb], 1, c, s);
                        a[((jj + 1) - 1) + (j - 1) * lda] = COMPLEX(c);
                        b[((jj + 1) - 1) + (j - 1) * ldb] = -conj(s);
                    }
                }
                //
                //              Update A by transformations from right.
                //
                jj = mod(ihi - j - 1, 3);
                for (i = ihi - j - 3; i >= jj + 1; i = i - 3) {
                    ctemp = a[((j + 1 + i) - 1) + (j - 1) * lda];
                    s = -b[((j + 1 + i) - 1) + (j - 1) * ldb];
                    c1 = a[((j + 2 + i) - 1) + (j - 1) * lda];
                    s1 = -b[((j + 2 + i) - 1) + (j - 1) * ldb];
                    c2 = a[((j + 3 + i) - 1) + (j - 1) * lda];
                    s2 = -b[((j + 3 + i) - 1) + (j - 1) * ldb];
                    //
                    for (k = top + 1; k <= ihi; k = k + 1) {
                        temp = a[(k - 1) + ((j + i) - 1) * lda];
                        temp1 = a[(k - 1) + ((j + i + 1) - 1) * lda];
                        temp2 = a[(k - 1) + ((j + i + 2) - 1) * lda];
                        temp3 = a[(k - 1) + ((j + i + 3) - 1) * lda];
                        a[(k - 1) + ((j + i + 3) - 1) * lda] = c2 * temp3 + conj(s2) * temp2;
                        temp2 = -s2 * temp3 + c2 * temp2;
                        a[(k - 1) + ((j + i + 2) - 1) * lda] = c1 * temp2 + conj(s1) * temp1;
                        temp1 = -s1 * temp2 + c1 * temp1;
                        a[(k - 1) + ((j + i + 1) - 1) * lda] = ctemp * temp1 + conj(s) * temp;
                        a[(k - 1) + ((j + i) - 1) * lda] = -s * temp1 + ctemp * temp;
                    }
                }
                //
                if (jj > 0) {
                    for (i = jj; i >= 1; i = i - 1) {
                        c = a[((j + 1 + i) - 1) + (j - 1) * lda].real();
                        Crot(ihi - top, a[((top + 1) - 1) + ((j + i + 1) - 1) * lda], 1, a[((top + 1) - 1) + ((j + i) - 1) * lda], 1, c, -conj(b[((j + 1 + i) - 1) + (j - 1) * ldb]));
                    }
                }
                //
                //              Update (J+1)th column of A by transformations from left.
                //
                if (j < jcol + nnb - 1) {
                    len = 1 + j - jcol;
                    //
                    //                 Multiply with the trailing accumulated unitary
                    //                 matrix, which takes the form
                    //
                    //                        [  U11  U12  ]
                    //                    U = [            ],
                    //                        [  U21  U22  ]
                    //
                    //                 where U21 is a LEN-by-LEN matrix and U12 is lower
                    //                 triangular.
                    //
                    jrow = ihi - nblst + 1;
                    Cgemv("Conjugate", nblst, len, cone, work, nblst, a[(jrow - 1) + ((j + 1) - 1) * lda], 1, czero, work[pw - 1], 1);
                    ppw = pw + len;
                    for (i = jrow; i <= jrow + nblst - len - 1; i = i + 1) {
                        work[ppw - 1] = a[(i - 1) + ((j + 1) - 1) * lda];
                        ppw++;
                    }
                    Ctrmv("Lower", "Conjugate", "Non-unit", nblst - len, work[(len * nblst + 1) - 1], nblst, work[(pw + len) - 1], 1);
                    Cgemv("Conjugate", len, nblst - len, cone, work[((len + 1) * nblst - len + 1) - 1], nblst, a[((jrow + nblst - len) - 1) + ((j + 1) - 1) * lda], 1, cone, work[(pw + len) - 1], 1);
                    ppw = pw;
                    for (i = jrow; i <= jrow + nblst - 1; i = i + 1) {
                        a[(i - 1) + ((j + 1) - 1) * lda] = work[ppw - 1];
                        ppw++;
                    }
                    //
                    //                 Multiply with the other accumulated unitary
                    //                 matrices, which take the form
                    //
                    //                        [  U11  U12   0  ]
                    //                        [                ]
                    //                    U = [  U21  U22   0  ],
                    //                        [                ]
                    //                        [   0    0    I  ]
                    //
                    //                 where I denotes the (NNB-LEN)-by-(NNB-LEN) identity
                    //                 matrix, U21 is a LEN-by-LEN upper triangular matrix
                    //                 and U12 is an NNB-by-NNB lower triangular matrix.
                    //
                    ppwo = 1 + nblst * nblst;
                    j0 = jrow - nnb;
                    for (jrow = j0; jrow <= jcol + 1; jrow = jrow + -nnb) {
                        ppw = pw + len;
                        for (i = jrow; i <= jrow + nnb - 1; i = i + 1) {
                            work[ppw - 1] = a[(i - 1) + ((j + 1) - 1) * lda];
                            ppw++;
                        }
                        ppw = pw;
                        for (i = jrow + nnb; i <= jrow + nnb + len - 1; i = i + 1) {
                            work[ppw - 1] = a[(i - 1) + ((j + 1) - 1) * lda];
                            ppw++;
                        }
                        Ctrmv("Upper", "Conjugate", "Non-unit", len, work[(ppwo + nnb) - 1], 2 * nnb, work[pw - 1], 1);
                        Ctrmv("Lower", "Conjugate", "Non-unit", nnb, work[(ppwo + 2 * len * nnb) - 1], 2 * nnb, work[(pw + len) - 1], 1);
                        Cgemv("Conjugate", nnb, len, cone, work[ppwo - 1], 2 * nnb, a[(jrow - 1) + ((j + 1) - 1) * lda], 1, cone, work[pw - 1], 1);
                        Cgemv("Conjugate", len, nnb, cone, work[(ppwo + 2 * len * nnb + nnb) - 1], 2 * nnb, a[((jrow + nnb) - 1) + ((j + 1) - 1) * lda], 1, cone, work[(pw + len) - 1], 1);
                        ppw = pw;
                        for (i = jrow; i <= jrow + len + nnb - 1; i = i + 1) {
                            a[(i - 1) + ((j + 1) - 1) * lda] = work[ppw - 1];
                            ppw++;
                        }
                        ppwo += 4 * nnb * nnb;
                    }
                }
            }
            //
            //           Apply accumulated unitary matrices to A.
            //
            cola = n - jcol - nnb + 1;
            j = ihi - nblst + 1;
            Cgemm("Conjugate", "No Transpose", nblst, cola, nblst, cone, work, nblst, a[(j - 1) + ((jcol + nnb) - 1) * lda], lda, czero, work[pw - 1], nblst);
            Clacpy("All", nblst, cola, work[pw - 1], nblst, a[(j - 1) + ((jcol + nnb) - 1) * lda], lda);
            ppwo = nblst * nblst + 1;
            j0 = j - nnb;
            for (j = j0; j <= jcol + 1; j = j + -nnb) {
                if (blk22) {
                    //
                    //                 Exploit the structure of
                    //
                    //                        [  U11  U12  ]
                    //                    U = [            ]
                    //                        [  U21  U22  ],
                    //
                    //                 where all blocks are NNB-by-NNB, U21 is upper
                    //                 triangular and U12 is lower triangular.
                    //
                    Cunm22("Left", "Conjugate", 2 * nnb, cola, nnb, nnb, work[ppwo - 1], 2 * nnb, a[(j - 1) + ((jcol + nnb) - 1) * lda], lda, work[pw - 1], lwork - pw + 1, ierr);
                } else {
                    //
                    //                 Ignore the structure of U.
                    //
                    Cgemm("Conjugate", "No Transpose", 2 * nnb, cola, 2 * nnb, cone, work[ppwo - 1], 2 * nnb, a[(j - 1) + ((jcol + nnb) - 1) * lda], lda, czero, work[pw - 1], 2 * nnb);
                    Clacpy("All", 2 * nnb, cola, work[pw - 1], 2 * nnb, a[(j - 1) + ((jcol + nnb) - 1) * lda], lda);
                }
                ppwo += 4 * nnb * nnb;
            }
            //
            //           Apply accumulated unitary matrices to Q.
            //
            if (wantq) {
                j = ihi - nblst + 1;
                if (initq) {
                    topq = max(2, j - jcol + 1);
                    nh = ihi - topq + 1;
                } else {
                    topq = 1;
                    nh = n;
                }
                Cgemm("No Transpose", "No Transpose", nh, nblst, nblst, cone, q[(topq - 1) + (j - 1) * ldq], ldq, work, nblst, czero, work[pw - 1], nh);
                Clacpy("All", nh, nblst, work[pw - 1], nh, q[(topq - 1) + (j - 1) * ldq], ldq);
                ppwo = nblst * nblst + 1;
                j0 = j - nnb;
                for (j = j0; j <= jcol + 1; j = j + -nnb) {
                    if (initq) {
                        topq = max(2, j - jcol + 1);
                        nh = ihi - topq + 1;
                    }
                    if (blk22) {
                        //
                        //                    Exploit the structure of U.
                        //
                        Cunm22("Right", "No Transpose", nh, 2 * nnb, nnb, nnb, work[ppwo - 1], 2 * nnb, q[(topq - 1) + (j - 1) * ldq], ldq, work[pw - 1], lwork - pw + 1, ierr);
                    } else {
                        //
                        //                    Ignore the structure of U.
                        //
                        Cgemm("No Transpose", "No Transpose", nh, 2 * nnb, 2 * nnb, cone, q[(topq - 1) + (j - 1) * ldq], ldq, work[ppwo - 1], 2 * nnb, czero, work[pw - 1], nh);
                        Clacpy("All", nh, 2 * nnb, work[pw - 1], nh, q[(topq - 1) + (j - 1) * ldq], ldq);
                    }
                    ppwo += 4 * nnb * nnb;
                }
            }
            //
            //           Accumulate right Givens rotations if required.
            //
            if (wantz || top > 0) {
                //
                //              Initialize small unitary factors that will hold the
                //              accumulated Givens rotations in workspace.
                //
                Claset("All", nblst, nblst, czero, cone, work, nblst);
                pw = nblst * nblst + 1;
                for (i = 1; i <= n2nb; i = i + 1) {
                    Claset("All", 2 * nnb, 2 * nnb, czero, cone, work[pw - 1], 2 * nnb);
                    pw += 4 * nnb * nnb;
                }
                //
                //              Accumulate Givens rotations INTEGERo workspace array.
                //
                for (j = jcol; j <= jcol + nnb - 1; j = j + 1) {
                    ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
                    len = 2 + j - jcol;
                    jrow = j + n2nb * nnb + 2;
                    for (i = ihi; i >= jrow; i = i - 1) {
                        ctemp = a[(i - 1) + (j - 1) * lda];
                        a[(i - 1) + (j - 1) * lda] = czero;
                        s = b[(i - 1) + (j - 1) * ldb];
                        b[(i - 1) + (j - 1) * ldb] = czero;
                        for (jj = ppw; jj <= ppw + len - 1; jj = jj + 1) {
                            temp = work[(jj + nblst) - 1];
                            work[(jj + nblst) - 1] = ctemp * temp - conj(s) * work[jj - 1];
                            work[jj - 1] = s * temp + ctemp * work[jj - 1];
                        }
                        len++;
                        ppw = ppw - nblst - 1;
                    }
                    //
                    ppwo = nblst * nblst + (nnb + j - jcol - 1) * 2 * nnb + nnb;
                    j0 = jrow - nnb;
                    for (jrow = j0; jrow <= j + 2; jrow = jrow + -nnb) {
                        ppw = ppwo;
                        len = 2 + j - jcol;
                        for (i = jrow + nnb - 1; i >= jrow; i = i - 1) {
                            ctemp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = czero;
                            s = b[(i - 1) + (j - 1) * ldb];
                            b[(i - 1) + (j - 1) * ldb] = czero;
                            for (jj = ppw; jj <= ppw + len - 1; jj = jj + 1) {
                                temp = work[(jj + 2 * nnb) - 1];
                                work[(jj + 2 * nnb) - 1] = ctemp * temp - conj(s) * work[jj - 1];
                                work[jj - 1] = s * temp + ctemp * work[jj - 1];
                            }
                            len++;
                            ppw = ppw - 2 * nnb - 1;
                        }
                        ppwo += 4 * nnb * nnb;
                    }
                }
            } else {
                //
                Claset("Lower", ihi - jcol - 1, nnb, czero, czero, a[((jcol + 2) - 1) + (jcol - 1) * lda], lda);
                Claset("Lower", ihi - jcol - 1, nnb, czero, czero, b[((jcol + 2) - 1) + (jcol - 1) * ldb], ldb);
            }
            //
            //           Apply accumulated unitary matrices to A and B.
            //
            if (top > 0) {
                j = ihi - nblst + 1;
                Cgemm("No Transpose", "No Transpose", top, nblst, nblst, cone, a[(j - 1) * lda], lda, work, nblst, czero, work[pw - 1], top);
                Clacpy("All", top, nblst, work[pw - 1], top, a[(j - 1) * lda], lda);
                ppwo = nblst * nblst + 1;
                j0 = j - nnb;
                for (j = j0; j <= jcol + 1; j = j + -nnb) {
                    if (blk22) {
                        //
                        //                    Exploit the structure of U.
                        //
                        Cunm22("Right", "No Transpose", top, 2 * nnb, nnb, nnb, work[ppwo - 1], 2 * nnb, a[(j - 1) * lda], lda, work[pw - 1], lwork - pw + 1, ierr);
                    } else {
                        //
                        //                    Ignore the structure of U.
                        //
                        Cgemm("No Transpose", "No Transpose", top, 2 * nnb, 2 * nnb, cone, a[(j - 1) * lda], lda, work[ppwo - 1], 2 * nnb, czero, work[pw - 1], top);
                        Clacpy("All", top, 2 * nnb, work[pw - 1], top, a[(j - 1) * lda], lda);
                    }
                    ppwo += 4 * nnb * nnb;
                }
                //
                j = ihi - nblst + 1;
                Cgemm("No Transpose", "No Transpose", top, nblst, nblst, cone, b[(j - 1) * ldb], ldb, work, nblst, czero, work[pw - 1], top);
                Clacpy("All", top, nblst, work[pw - 1], top, b[(j - 1) * ldb], ldb);
                ppwo = nblst * nblst + 1;
                j0 = j - nnb;
                for (j = j0; j <= jcol + 1; j = j + -nnb) {
                    if (blk22) {
                        //
                        //                    Exploit the structure of U.
                        //
                        Cunm22("Right", "No Transpose", top, 2 * nnb, nnb, nnb, work[ppwo - 1], 2 * nnb, b[(j - 1) * ldb], ldb, work[pw - 1], lwork - pw + 1, ierr);
                    } else {
                        //
                        //                    Ignore the structure of U.
                        //
                        Cgemm("No Transpose", "No Transpose", top, 2 * nnb, 2 * nnb, cone, b[(j - 1) * ldb], ldb, work[ppwo - 1], 2 * nnb, czero, work[pw - 1], top);
                        Clacpy("All", top, 2 * nnb, work[pw - 1], top, b[(j - 1) * ldb], ldb);
                    }
                    ppwo += 4 * nnb * nnb;
                }
            }
            //
            //           Apply accumulated unitary matrices to Z.
            //
            if (wantz) {
                j = ihi - nblst + 1;
                if (initq) {
                    topq = max(2, j - jcol + 1);
                    nh = ihi - topq + 1;
                } else {
                    topq = 1;
                    nh = n;
                }
                Cgemm("No Transpose", "No Transpose", nh, nblst, nblst, cone, z[(topq - 1) + (j - 1) * ldz], ldz, work, nblst, czero, work[pw - 1], nh);
                Clacpy("All", nh, nblst, work[pw - 1], nh, z[(topq - 1) + (j - 1) * ldz], ldz);
                ppwo = nblst * nblst + 1;
                j0 = j - nnb;
                for (j = j0; j <= jcol + 1; j = j + -nnb) {
                    if (initq) {
                        topq = max(2, j - jcol + 1);
                        nh = ihi - topq + 1;
                    }
                    if (blk22) {
                        //
                        //                    Exploit the structure of U.
                        //
                        Cunm22("Right", "No Transpose", nh, 2 * nnb, nnb, nnb, work[ppwo - 1], 2 * nnb, z[(topq - 1) + (j - 1) * ldz], ldz, work[pw - 1], lwork - pw + 1, ierr);
                    } else {
                        //
                        //                    Ignore the structure of U.
                        //
                        Cgemm("No Transpose", "No Transpose", nh, 2 * nnb, 2 * nnb, cone, z[(topq - 1) + (j - 1) * ldz], ldz, work[ppwo - 1], 2 * nnb, czero, work[pw - 1], nh);
                        Clacpy("All", nh, 2 * nnb, work[pw - 1], nh, z[(topq - 1) + (j - 1) * ldz], ldz);
                    }
                    ppwo += 4 * nnb * nnb;
                }
            }
        }
    }
    //
    //     Use unblocked code to reduce the rest of the matrix
    //     Avoid re-initialization of modified Q and Z.
    //
    str<1> compq2 = compq;
    str<1> compz2 = compz;
    if (jcol != ilo) {
        if (wantq) {
            compq2 = "V";
        }
        if (wantz) {
            compz2 = "V";
        }
    }
    //
    if (jcol < ihi) {
        Cgghrd(compq2, compz2, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, ierr);
    }
    work[1 - 1] = COMPLEX(lwkopt);
    //
    //     End of Cgghd3
    //
}
