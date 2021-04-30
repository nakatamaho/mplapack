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

INTEGER iMparam2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const ni, INTEGER const nbi, INTEGER const ibi, INTEGER const nxi) {
    INTEGER return_value = 0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //
    //  ================================================================
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Invalid value for ISPEC
    //
    if ((ispec < 17) || (ispec > 21)) {
        return_value = -1;
        return return_value;
    }
    //
    //     Get the number of threads
    //
    INTEGER nthreads = 1;
    //      WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC
    //
    char subnam[12];
    INTEGER ic = 0;
    INTEGER iz = 0;
    INTEGER i = 0;
    char prec;
    char algo[3];
    char stag[3];
    bool rprec = false;
    bool cprec = false;
    if (ispec != 19) {
        //
        //        Convert NAME to upper case if the first character is lower case.
        //
        return_value = -1;
        subnam = name;
        ic = ichar(subnam[(1 - 1)]);
        iz = ichar("Z");
        if (iz == 90 || iz == 122) {
            //
            //           ASCII character set
            //
            if (ic >= 97 && ic <= 122) {
                subnam[(1 - 1)] = fchar[(ic - 32) - 1];
                for (i = 2; i <= 12; i = i + 1) {
                    ic = ichar(subnam[(i - 1) + (i - 1) * ldsubnam]);
                    if (ic >= 97 && ic <= 122) {
                        subnam[(i - 1) + (i - 1) * ldsubnam] = fchar[(ic - 32) - 1];
                    }
                }
            }
            //
        } else if (iz == 233 || iz == 169) {
            //
            //           EBCDIC character set
            //
            if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
                subnam[(1 - 1)] = fchar[(ic + 64) - 1];
                for (i = 2; i <= 12; i = i + 1) {
                    ic = ichar(subnam[(i - 1) + (i - 1) * ldsubnam]);
                    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
                        subnam[(i - 1) + (i - 1) * ldsubnam] = fchar[(ic + 64) - 1];
                    }
                }
            }
            //
        } else if (iz == 218 || iz == 250) {
            //
            //           Prime machines:  ASCII+128
            //
            if (ic >= 225 && ic <= 250) {
                subnam[(1 - 1)] = fchar[(ic - 32) - 1];
                for (i = 2; i <= 12; i = i + 1) {
                    ic = ichar(subnam[(i - 1) + (i - 1) * ldsubnam]);
                    if (ic >= 225 && ic <= 250) {
                        subnam[(i - 1) + (i - 1) * ldsubnam] = fchar[(ic - 32) - 1];
                    }
                }
            }
        }
        //
        prec = subnam[(1 - 1)];
        algo = subnam[(4 - 1) + (6 - 1) * ldsubnam];
        stag = subnam[(8 - 1) + (12 - 1) * ldsubnam];
        rprec = prec == "R";
        cprec = prec == "C";
        //
        //        Invalid value for PRECISION
        //
        if (!(rprec || cprec)) {
            return_value = -1;
            return return_value;
        }
    }
    //      WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC,
    //     $           '   ALGO ',ALGO,'    STAGE ',STAG
    //
    INTEGER kd = 0;
    INTEGER ib = 0;
    char vect;
    INTEGER lhous = 0;
    INTEGER lwork = 0;
    INTEGER qroptnb = 0;
    INTEGER lqoptnb = 0;
    INTEGER factoptnb = 0;
    if ((ispec == 17) || (ispec == 18)) {
        //
        //     ISPEC = 17, 18:  block size KD, IB
        //     Could be also dependent from N but for now it
        //     depend only on sequential or parallel
        //
        if (nthreads > 4) {
            if (cprec) {
                kd = 128;
                ib = 32;
            } else {
                kd = 160;
                ib = 40;
            }
        } else if (nthreads > 1) {
            if (cprec) {
                kd = 64;
                ib = 32;
            } else {
                kd = 64;
                ib = 32;
            }
        } else {
            if (cprec) {
                kd = 16;
                ib = 16;
            } else {
                kd = 32;
                ib = 16;
            }
        }
        if (ispec == 17) {
            return_value = kd;
        }
        if (ispec == 18) {
            return_value = ib;
        }
        //
    } else if (ispec == 19) {
        //
        //     ISPEC = 19:
        //     LHOUS length of the Houselholder representation
        //     matrix (V,T) of the second stage. should be >= 1.
        //
        //     Will add the VECT OPTION HERE next release
        vect = opts[(1 - 1)];
        if (vect == "N") {
            lhous = max((INTEGER)1, 4 * ni);
        } else {
            //           This is not correct, it need to call the ALGO and the stage2
            lhous = max((INTEGER)1, 4 * ni) + ibi;
        }
        if (lhous >= 0) {
            return_value = lhous;
        } else {
            return_value = -1;
        }
        //
    } else if (ispec == 20) {
        //
        //     ISPEC = 20: (21 for future use)
        //     LWORK length of the workspace for
        //     either or both stages for TRD and BRD. should be >= 1.
        //     TRD:
        //     TRD_stage 1: = LT + LW + LS1 + LS2
        //                  = LDT*KD + N*KD + N*MAX(KD,FACTOPTNB) + LDS2*KD
        //                    where LDT=LDS2=KD
        //                  = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
        //     TRD_stage 2: = (2NB+1)*N + KD*NTHREADS
        //     TRD_both   : = max(stage1,stage2) + AB ( AB=(KD+1)*N )
        //                  = N*KD + N*max(KD+1,FACTOPTNB)
        //                    + max((INTEGER)2*KD*KD, KD*NTHREADS)
        //                    + (KD+1)*N
        lwork = -1;
        subnam[(1 - 1)] = prec;
        subnam[(2 - 1) + (6 - 1) * ldsubnam] = "GEQRF";
        qroptnb = iMlaenv(1, subnam, " ", ni, nbi, -1, -1);
        subnam[(2 - 1) + (6 - 1) * ldsubnam] = "GELQF";
        lqoptnb = iMlaenv(1, subnam, " ", nbi, ni, -1, -1);
        //        Could be QR or LQ for TRD and the max for BRD
        factoptnb = max(qroptnb, lqoptnb);
        if (algo == "TRD") {
            if (stag == "2STAG") {
                lwork = ni * nbi + ni * max(nbi + 1, factoptnb) + max((INTEGER)2 * nbi * nbi, nbi * nthreads) + (nbi + 1) * ni;
            } else if ((stag == "HE2HB") || (stag == "SY2SB")) {
                lwork = ni * nbi + ni * max(nbi, factoptnb) + 2 * nbi * nbi;
            } else if ((stag == "HB2ST") || (stag == "SB2ST")) {
                lwork = (2 * nbi + 1) * ni + nbi * nthreads;
            }
        } else if (algo == "BRD") {
            if (stag == "2STAG") {
                lwork = 2 * ni * nbi + ni * max(nbi + 1, factoptnb) + max((INTEGER)2 * nbi * nbi, nbi * nthreads) + (nbi + 1) * ni;
            } else if (stag == "GE2GB") {
                lwork = ni * nbi + ni * max(nbi, factoptnb) + 2 * nbi * nbi;
            } else if (stag == "GB2BD") {
                lwork = (3 * nbi + 1) * ni + nbi * nthreads;
            }
        }
        lwork = max((INTEGER)1, lwork);
        //
        if (lwork > 0) {
            return_value = lwork;
        } else {
            return_value = -1;
        }
        //
    } else if (ispec == 21) {
        //
        //     ISPEC = 21 for future use
        return_value = nxi;
    }
    return return_value;
    //
    //     ==== End of iMparam2stage ====
    //
}
