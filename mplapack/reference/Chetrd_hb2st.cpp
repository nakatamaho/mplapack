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

void Chetrd_hb2st(const char *stage1, const char *vect, const char *uplo, INTEGER const n, INTEGER const kd, COMPLEX *ab, INTEGER const ldab, REAL *d, REAL *e, COMPLEX *hous, INTEGER const lhous, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //     Determine the minimal workspace size required.
    //     Test the input parameters
    //
    INTEGER debug = 0;
    info = 0;
    bool afters1 = Mlsame(stage1, "Y");
    bool wantq = Mlsame(vect, "V");
    bool upper = Mlsame(uplo, "U");
    bool lquery = (lwork == -1) || (lhous == -1);
    //
    //     Determine the block size, the workspace size and the hous size.
    //
    INTEGER ib = iMlaenv2stage(2, "Chetrd_hb2st", vect, n, kd, -1, -1);
    INTEGER lhmin = iMlaenv2stage(3, "Chetrd_hb2st", vect, n, kd, ib, -1);
    INTEGER lwmin = iMlaenv2stage(4, "Chetrd_hb2st", vect, n, kd, ib, -1);
    //
    if (!afters1 && !Mlsame(stage1, "N")) {
        info = -1;
    } else if (!Mlsame(vect, "N")) {
        info = -2;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (kd < 0) {
        info = -5;
    } else if (ldab < (kd + 1)) {
        info = -7;
    } else if (lhous < lhmin && !lquery) {
        info = -11;
    } else if (lwork < lwmin && !lquery) {
        info = -13;
    }
    //
    if (info == 0) {
        hous[1 - 1] = lhmin;
        work[1 - 1] = lwmin;
    }
    //
    if (info != 0) {
        Mxerbla("Chetrd_hb2st", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        hous[1 - 1] = 1;
        work[1 - 1] = 1;
        return;
    }
    //
    //     Determine pointer position
    //
    INTEGER ldv = kd + ib;
    INTEGER sizetau = 2 * n;
    INTEGER sizev = 2 * n;
    INTEGER indtau = 1;
    INTEGER indv = indtau + sizetau;
    INTEGER lda = 2 * kd + 1;
    INTEGER sizea = lda * n;
    INTEGER inda = 1;
    INTEGER indw = inda + sizea;
    INTEGER nthreads = 1;
    INTEGER tid = 0;
    //
    INTEGER apos = 0;
    INTEGER awpos = 0;
    INTEGER dpos = 0;
    INTEGER ofdpos = 0;
    INTEGER abdpos = 0;
    INTEGER abofdpos = 0;
    if (upper) {
        apos = inda + kd;
        awpos = inda;
        dpos = apos + kd;
        ofdpos = dpos - 1;
        abdpos = kd + 1;
        abofdpos = kd;
    } else {
        apos = inda;
        awpos = inda + kd + 1;
        dpos = apos;
        ofdpos = dpos + 1;
        abdpos = 1;
        abofdpos = 2;
        //
    }
    //
    //     Case KD=0:
    //     The matrix is diagonal. We just copy it (convert to "real" for
    //     complex because D is REAL and the imaginary part should be 0)
    //     and store it in D. A sequential code here is better or
    //     in a parallel environment it might need two cores for D and E
    //
    INTEGER i = 0;
    const REAL rzero = 0.0;
    if (kd == 0) {
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = ab[(abdpos - 1) + (i - 1) * ldab].real();
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            e[i - 1] = rzero;
        }
        //
        hous[1 - 1] = 1;
        work[1 - 1] = 1;
        return;
    }
    //
    //     Case KD=1:
    //     The matrix is already Tridiagonal. We have to make diagonal
    //     and offdiagonal elements real, and store them in D and E.
    //     For that, for real precision just copy the diag and offdiag
    //     to D and E while for the COMPLEX case the bulge chasing is
    //     performed to convert the hermetian tridiagonal to symmetric
    //     tridiagonal. A simpler conversion formula might be used, but then
    //     updating the Q matrix will be required and based if Q is generated
    //     or not this might complicate the story.
    //
    COMPLEX tmp = 0.0;
    REAL abstmp = 0.0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    if (kd == 1) {
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = ab[(abdpos - 1) + (i - 1) * ldab].real();
        }
        //
        //         make off-diagonal elements real and copy them to E
        //
        if (upper) {
            for (i = 1; i <= n - 1; i = i + 1) {
                tmp = ab[(abofdpos - 1) + ((i + 1) - 1) * ldab];
                abstmp = abs(tmp);
                ab[(abofdpos - 1) + ((i + 1) - 1) * ldab] = abstmp;
                e[i - 1] = abstmp;
                if (abstmp != rzero) {
                    tmp = tmp / abstmp;
                } else {
                    tmp = one;
                }
                if (i < n - 1) {
                    ab[(abofdpos - 1) + ((i + 2) - 1) * ldab] = ab[(abofdpos - 1) + ((i + 2) - 1) * ldab] * tmp;
                }
                //                  IF( WANTZ ) THEN
                //                     CALL Cscal( N, DCONJG( TMP ), Q( 1, I+1 ), 1 )
                //                  END IF
            }
        } else {
            for (i = 1; i <= n - 1; i = i + 1) {
                tmp = ab[(abofdpos - 1) + (i - 1) * ldab];
                abstmp = abs(tmp);
                ab[(abofdpos - 1) + (i - 1) * ldab] = abstmp;
                e[i - 1] = abstmp;
                if (abstmp != rzero) {
                    tmp = tmp / abstmp;
                } else {
                    tmp = one;
                }
                if (i < n - 1) {
                    ab[(abofdpos - 1) + ((i + 1) - 1) * ldab] = ab[(abofdpos - 1) + ((i + 1) - 1) * ldab] * tmp;
                }
                //                 IF( WANTQ ) THEN
                //                    CALL Cscal( N, TMP, Q( 1, I+1 ), 1 )
                //                 END IF
            }
        }
        //
        hous[1 - 1] = 1;
        work[1 - 1] = 1;
        return;
    }
    //
    //     Main code start here.
    //     Reduce the hermitian band of A to a tridiagonal matrix.
    //
    INTEGER thgrsiz = n;
    INTEGER grsiz = 1;
    INTEGER shift = 3;
    INTEGER nbtiles = ceil(double(n) / double(kd));
    INTEGER stepercol = ceil(double(shift) / double(grsiz));
    INTEGER thgrnb = ceil(double(n - 1) / double(thgrsiz));
    //
    Clacpy("A", kd + 1, n, ab, ldab, &work[apos - 1], lda);
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    Claset("A", kd, n, zero, zero, &work[awpos - 1], lda);
    //
    //     main bulge chasing loop
    //
    INTEGER thgrid = 0;
    INTEGER stt = 0;
    INTEGER thed = 0;
    INTEGER ed = 0;
    INTEGER m = 0;
    INTEGER st = 0;
    INTEGER sweepid = 0;
    INTEGER k = 0;
    INTEGER myid = 0;
    INTEGER ttype = 0;
    INTEGER colpt = 0;
    INTEGER stind = 0;
    INTEGER edind = 0;
    INTEGER blklastind = 0;
    for (thgrid = 1; thgrid <= thgrnb; thgrid = thgrid + 1) {
        stt = (thgrid - 1) * thgrsiz + 1;
        thed = min((stt + thgrsiz - 1), (n - 1));
        for (i = stt; i <= n - 1; i = i + 1) {
            ed = min(i, thed);
            if (stt > ed) {
                break;
            }
            for (m = 1; m <= stepercol; m = m + 1) {
                st = stt;
                for (sweepid = st; sweepid <= ed; sweepid = sweepid + 1) {
                    for (k = 1; k <= grsiz; k = k + 1) {
                        myid = (i - sweepid) * (stepercol * grsiz) + (m - 1) * grsiz + k;
                        if (myid == 1) {
                            ttype = 1;
                        } else {
                            ttype = mod(myid, 2) + 2;
                        }
                        //
                        if (ttype == 2) {
                            colpt = (myid / 2) * kd + sweepid;
                            stind = colpt - kd + 1;
                            edind = min(colpt, n);
                            blklastind = colpt;
                        } else {
                            colpt = ((myid + 1) / 2) * kd + sweepid;
                            stind = colpt - kd + 1;
                            edind = min(colpt, n);
                            if ((stind >= edind - 1) && (edind == n)) {
                                blklastind = n;
                            } else {
                                blklastind = 0;
                            }
                        }
                        //
                        //                         Call the kernel
                        //
                        Chb2st_kernels(uplo, wantq, ttype, stind, edind, sweepid, n, kd, ib, &work[inda - 1], lda, &hous[indv - 1], &hous[indtau - 1], ldv, &work[(indw + tid * kd) - 1]);
                        if (blklastind >= (n - 1)) {
                            stt++;
                            break;
                        }
                    }
                }
            }
        }
    }
    //
    //     Copy the diagonal from A to D. Note that D is REAL thus only
    //     the Real part is needed, the imaginary part should be zero.
    //
    for (i = 1; i <= n; i = i + 1) {
        d[i - 1] = work[(dpos + (i - 1) * lda) - 1].real();
    }
    //
    //     Copy the off diagonal from A to E. Note that E is REAL thus only
    //     the Real part is needed, the imaginary part should be zero.
    //
    if (upper) {
        for (i = 1; i <= n - 1; i = i + 1) {
            e[i - 1] = work[(ofdpos + i * lda) - 1].real();
        }
    } else {
        for (i = 1; i <= n - 1; i = i + 1) {
            e[i - 1] = work[(ofdpos + (i - 1) * lda) - 1].real();
        }
    }
    //
    hous[1 - 1] = lhmin;
    work[1 - 1] = lwmin;
    //
    //     End of Chetrd_hb2st
    //
}
