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

void Alaerh(common &cmn, const char *path, const char *subnam, INTEGER const info, INTEGER const infoe, const char *opts, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const n5, INTEGER const imat, INTEGER const nfail, INTEGER &nerrs, INTEGER const nout) {
    common_write write(cmn);
    static const char *format_9949 = "(' ==> Doing only the condition estimate for this case')";
    static const char *format_9952 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', TRANS=''',a1,"
                                     "''', DIAG=''',a1,''', NORMIN=''',a1,''', N =',i5,', type ',i2)";
    static const char *format_9953 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', TRANS=''',a1,"
                                     "''', DIAG=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    static const char *format_9955 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', NRHS =',i4,', type ',i2)";
    static const char *format_9956 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', NB =',i4,', type ',i2)";
    static const char *format_9957 = "(' *** Error code from ',a,'=',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', KD =',i5,', NRHS =',i4,', type ',i2)";
    static const char *format_9960 = "(' *** Error code from ',a,' =',i5,' for UPLO = ''',a1,''', N =',i5,"
                                     "', type ',i2)";
    static const char *format_9963 = "(' *** Error code from ',a,' =',i5,/,' ==> TRANS = ''',a1,''', N =',i5,"
                                     "', NRHS =',i4,', type ',i2)";
    static const char *format_9967 = "(' *** Error code from ',a,' =',i5,/,' ==> NORM=''',a1,''', UPLO =''',a1,"
                                     "''', DIAG=''',a1,''', N =',i5,', type ',i2)";
    static const char *format_9969 = "(' *** Error code from ',a,' =',i5,' for NORM = ''',a1,''', N =',i5,"
                                     "', type ',i2)";
    static const char *format_9970 = "(' *** Error code from ',a,' =',i5,' for N =',i5,', NRHS =',i4,', type ',"
                                     "i2)";
    static const char *format_9971 = "(' *** Error code from ',a,'=',i5,' for N=',i5,', NB=',i4,', type ',i2)";
    static const char *format_9973 = "(' *** Error code from ',a,' =',i5,' for N =',i5,', type ',i2)";
    static const char *format_9974 = "(' *** Error code from ',a,'=',i5,/,' ==> M =',i5,', N =',i5,', NRHS =',"
                                     "i4,', NB =',i4,', type ',i2)";
    static const char *format_9975 = "(' *** Error code from ',a,'=',i5,' for M=',i5,', N=',i5,', NB=',i4,"
                                     "', type ',i2)";
    static const char *format_9978 = "(' *** Error code from ',a,' =',i5,' for M =',i5,', N =',i5,', type ',i2)";
    static const char *format_9979 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                                     "' ==> UPLO = ''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    static const char *format_9980 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                                     "' ==> UPLO = ''',a1,''', N =',i5,', NB =',i4,', type ',i2)";
    static const char *format_9984 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> N =',i5,"
                                     "', NRHS =',i4,', type ',i2)";
    static const char *format_9987 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,' for N=',i5,"
                                     "', type ',i2)";
    static const char *format_9988 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> M =',i5,"
                                     "', N =',i5,', NB =',i4,', type ',i2)";
    static const char *format_9990 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', UPLO=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    static const char *format_9992 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', TRANS=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    static const char *format_9995 = "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,''', UPLO=''',a1,"
                                     "''', N =',i5,', NRHS =',i4,', type ',i2)";
    static const char *format_9997 = "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,''', TRANS=''',a1,"
                                     "''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    if (info == 0) {
        return;
    }
    str<2> p2 = path[(2 - 1) + (3 - 1) * ldpath];
    str<3> c3 = subnam[(4 - 1) + (6 - 1) * ldsubnam];
    //
    //     PrINTEGER the header if this is the first error message.
    //
    if (nfail == 0 && nerrs == 0) {
        if (Mlsamen3, c3, "SV " || Mlsamen3, c3, "SVX") {
            Aladhd(nout, path);
        } else {
            Alahd(nout, path);
        }
    }
    nerrs++;
    //
    //     PrINTEGER the message detailing the error and form of recovery,
    //     if any.
    //
    char uplo = char0;
    if (Mlsamen2, p2, "GE") {
        //
        //        xGE:  General matrices
        //
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, format_9988), subnam(1, len_trim(subnam)), info, infoe, m, n, n5, imat;
            } else {
                write(nout, format_9975), subnam(1, len_trim(subnam)), info, m, n, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9984), subnam(1, len_trim(subnam)), info, infoe, n, n5, imat;
            } else {
                write(nout, format_9970), subnam(1, len_trim(subnam)), info, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9992), subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9997), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "TRI") {
            //
            write(nout, format_9971), subnam(1, len_trim(subnam)), info, n, n5, imat;
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS") {
            //
            write(nout, format_9978), subnam(1, len_trim(subnam)), info, m, n, imat;
            //
        } else if (Mlsamen3, c3, "CON") {
            //
            write(nout, format_9969), subnam(1, len_trim(subnam)), info, opts(1, 1), m, imat;
            //
        } else if (Mlsamen3, c3, "LS ") {
            //
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> TRANS = ''',a1,''', M =',"
                        "i5,', N =',i5,', NRHS =',i4,', NB =',i4,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), m, n, kl, n5, imat;
            //
        } else if (Mlsamen3, c3, "LSX" || Mlsamen3, c3, "LSS") {
            //
            write(nout, format_9974), subnam(1, len_trim(subnam)), info, m, n, kl, n5, imat;
            //
        } else {
            //
            write(nout, format_9963), subnam(1, len_trim(subnam)), info, opts(1, 1), m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "GB") {
        //
        //        xGB:  General band matrices
        //
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> M = ',i5,', N =',i5,', KL =',i5,', KU =',i5,', NB =',i4,"
                            "', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, m, n, kl, ku, n5, imat;
            } else {
                write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> M = ',i5,', N =',i5,"
                            "', KL =',i5,', KU =',i5,', NB =',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, m, n, kl, ku, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> N =',i5,', KL =',i5,', KU =',i5,', NRHS =',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, n, kl, ku, n5, imat;
            } else {
                write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> N =',i5,', KL =',i5,"
                            "', KU =',i5,', NRHS =',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, n, kl, ku, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> FACT=''',a1,''', TRANS=''',a1,''', N=',i5,', KL=',i5,', KU=',"
                            "i5,', NRHS=',i4,', type ',i1)"),
                    subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, kl, ku, n5, imat;
            } else {
                write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,"
                            "''', TRANS=''',a1,''', N=',i5,', KL=',i5,', KU=',i5,', NRHS=',i4,"
                            "', type ',i1)"),
                    subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, kl, ku, n5, imat;
            }
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS") {
            //
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> M = ',i5,', N =',i5,"
                        "', KL =',i5,', KU =',i5,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, m, n, kl, ku, imat;
            //
        } else if (Mlsamen3, c3, "CON") {
            //
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> NORM =''',a1,''', N =',i5,"
                        "', KL =',i5,', KU =',i5,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), m, kl, ku, imat;
            //
        } else {
            //
            write(nout, "(' *** Error code from ',a,'=',i5,/,' ==> TRANS=''',a1,''', N =',i5,"
                        "', KL =',i5,', KU =',i5,', NRHS =',i4,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), m, kl, ku, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "GT") {
        //
        //        xGT:  General tridiagonal matrices
        //
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, format_9987), subnam(1, len_trim(subnam)), info, infoe, n, imat;
            } else {
                write(nout, format_9973), subnam(1, len_trim(subnam)), info, n, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9984), subnam(1, len_trim(subnam)), info, infoe, n, n5, imat;
            } else {
                write(nout, format_9970), subnam(1, len_trim(subnam)), info, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9992), subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9997), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "CON") {
            //
            write(nout, format_9969), subnam(1, len_trim(subnam)), info, opts(1, 1), m, imat;
            //
        } else {
            //
            write(nout, format_9963), subnam(1, len_trim(subnam)), info, opts(1, 1), m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "PO") {
        //
        //        xPO:  Symmetric or Hermitian positive definite matrices
        //
        uplo = opts[(1 - 1)];
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, format_9980), subnam(1, len_trim(subnam)), info, infoe, uplo, m, n5, imat;
            } else {
                write(nout, format_9956), subnam(1, len_trim(subnam)), info, uplo, m, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam(1, len_trim(subnam)), info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam(1, len_trim(subnam)), info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "TRI") {
            //
            write(nout, format_9956), subnam(1, len_trim(subnam)), info, uplo, m, n5, imat;
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS" || Mlsamen3, c3, "CON") {
            //
            write(nout, format_9960), subnam(1, len_trim(subnam)), info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam(1, len_trim(subnam)), info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "PS") {
        //
        //        xPS:  Symmetric or Hermitian positive semi-definite matrices
        //
        uplo = opts[(1 - 1)];
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, format_9980), subnam, info, infoe, uplo, m, n5, imat;
            } else {
                write(nout, format_9956), subnam, info, uplo, m, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam, info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam, info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam, info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam, info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "TRI") {
            //
            write(nout, format_9956), subnam, info, uplo, m, n5, imat;
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMT" || Mlsamen3, c3, "CON") {
            //
            write(nout, format_9960), subnam, info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam, info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "SY" || Mlsamen2, p2, "SR" || Mlsamen2, p2, "SK" || Mlsamen2, p2, "HE" || Mlsamen2, p2, "HR" || Mlsamen2, p2, "HK" || Mlsamen2, p2, "HA") {
        //
        //        xSY: symmetric indefinite matrices
        //             with partial (Bunch-Kaufman) pivoting;
        //        xSR: symmetric indefinite matrices
        //             with rook (bounded Bunch-Kaufman) pivoting;
        //        xSK: symmetric indefinite matrices
        //             with rook (bounded Bunch-Kaufman) pivoting,
        //             new storage format;
        //        xHE: Hermitian indefinite matrices
        //             with partial (Bunch-Kaufman) pivoting.
        //        xHR: Hermitian indefinite matrices
        //             with rook (bounded Bunch-Kaufman) pivoting;
        //        xHK: Hermitian indefinite matrices
        //             with rook (bounded Bunch-Kaufman) pivoting,
        //             new storage format;
        //        xHA: Hermitian matrices
        //             Aasen Algorithm
        //
        uplo = opts[(1 - 1)];
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, format_9980), subnam(1, len_trim(subnam)), info, infoe, uplo, m, n5, imat;
            } else {
                write(nout, format_9956), subnam(1, len_trim(subnam)), info, uplo, m, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen2, c3, "SV") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam(1, len_trim(subnam)), info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam(1, len_trim(subnam)), info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS" || Mlsamen3, c3, "TRI" || Mlsamen3, c3, "CON") {
            //
            write(nout, format_9960), subnam(1, len_trim(subnam)), info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam(1, len_trim(subnam)), info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "PP" || Mlsamen2, p2, "SP" || Mlsamen2, p2, "HP") {
        //
        //        xPP, xHP, or xSP:  Symmetric or Hermitian packed matrices
        //
        uplo = opts[(1 - 1)];
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> UPLO = ''',a1,''', N =',i5,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, uplo, m, imat;
            } else {
                write(nout, format_9960), subnam(1, len_trim(subnam)), info, uplo, m, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam(1, len_trim(subnam)), info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam(1, len_trim(subnam)), info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS" || Mlsamen3, c3, "TRI" || Mlsamen3, c3, "CON") {
            //
            write(nout, format_9960), subnam(1, len_trim(subnam)), info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam(1, len_trim(subnam)), info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "PB") {
        //
        //        xPB:  Symmetric (Hermitian) positive definite band matrix
        //
        uplo = opts[(1 - 1)];
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> UPLO = ''',a1,''', N =',i5,', KD =',i5,', NB =',i4,', type ',"
                            "i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, uplo, m, kl, n5, imat;
            } else {
                write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',"
                            "i5,', KD =',i5,', NB =',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, uplo, m, kl, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> UPLO=''',a1,''', N =',i5,', KD =',i5,', NRHS =',i4,', type ',"
                            "i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, uplo, n, kl, n5, imat;
            } else {
                write(nout, format_9957), subnam(1, len_trim(subnam)), info, uplo, n, kl, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> FACT=''',a1,''', UPLO=''',a1,''', N=',i5,', KD=',i5,"
                            "', NRHS=',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, kl, n5, imat;
            } else {
                write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,"
                            "''', UPLO=''',a1,''', N=',i5,', KD=',i5,', NRHS=',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, kl, n5, imat;
            }
            //
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS" || Mlsamen3, c3, "CON") {
            //
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',"
                        "i5,', KD =',i5,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, uplo, m, kl, imat;
            //
        } else {
            //
            write(nout, format_9957), subnam(1, len_trim(subnam)), info, uplo, m, kl, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "PT") {
        //
        //        xPT:  Positive definite tridiagonal matrices
        //
        if (Mlsamen3, c3, "TRF") {
            if (info != infoe && infoe != 0) {
                write(nout, format_9987), subnam(1, len_trim(subnam)), info, infoe, n, imat;
            } else {
                write(nout, format_9973), subnam(1, len_trim(subnam)), info, n, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen3, c3, "SV ") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9984), subnam(1, len_trim(subnam)), info, infoe, n, n5, imat;
            } else {
                write(nout, format_9970), subnam(1, len_trim(subnam)), info, n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "SVX") {
            //
            if (info != infoe && infoe != 0) {
                write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                            "' ==> FACT=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, infoe, opts(1, 1), n, n5, imat;
            } else {
                write(nout, "(' *** Error code from ',a,'=',i5,', FACT=''',a1,''', N=',i5,"
                            "', NRHS=',i4,', type ',i2)"),
                    subnam(1, len_trim(subnam)), info, opts(1, 1), n, n5, imat;
            }
            //
        } else if (Mlsamen3, c3, "CON") {
            //
            if (Mlsame(subnam[(1 - 1)], "S") || Mlsame(subnam[(1 - 1)], "D")) {
                write(nout, format_9973), subnam(1, len_trim(subnam)), info, m, imat;
            } else {
                write(nout, format_9969), subnam(1, len_trim(subnam)), info, opts(1, 1), m, imat;
            }
            //
        } else {
            //
            write(nout, format_9963), subnam(1, len_trim(subnam)), info, opts(1, 1), m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "TR") {
        //
        //        xTR:  Triangular matrix
        //
        if (Mlsamen3, c3, "TRI") {
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,"
                        "''', DIAG =''',a1,''', N =',i5,', NB =',i4,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), m, n5, imat;
        } else if (Mlsamen3, c3, "CON") {
            write(nout, format_9967), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATRS") {
            write(nout, format_9952), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), opts(4, 4), m, imat;
        } else {
            write(nout, format_9953), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "TP") {
        //
        //        xTP:  Triangular packed matrix
        //
        if (Mlsamen3, c3, "TRI") {
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,"
                        "''', DIAG =''',a1,''', N =',i5,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), m, imat;
        } else if (Mlsamen3, c3, "CON") {
            write(nout, format_9967), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATPS") {
            write(nout, format_9952), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), opts(4, 4), m, imat;
        } else {
            write(nout, format_9953), subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "TB") {
        //
        //        xTB:  Triangular band matrix
        //
        if (Mlsamen3, c3, "CON") {
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> NORM=''',a1,"
                        "''', UPLO =''',a1,''', DIAG=''',a1,''', N=',i5,', KD=',i5,', type ',"
                        "i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, kl, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATBS") {
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,"
                        "''', TRANS=''',a1,''', DIAG=''',a1,''', NORMIN=''',a1,''', N=',i5,"
                        "', KD=',i5,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), opts(4, 4), m, kl, imat;
        } else {
            write(nout, "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,"
                        "''', TRANS=''',a1,''', DIAG=''',a1,''', N=',i5,', KD=',i5,', NRHS=',"
                        "i4,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, kl, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "QR") {
        //
        //        xQR:  QR factorization
        //
        if (Mlsamen3, c3, "QRS") {
            write(nout, format_9974), subnam(1, len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS") {
            write(nout, format_9978), subnam(1, len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen2, p2, "LQ") {
        //
        //        xLQ:  LQ factorization
        //
        if (Mlsamen3, c3, "LQS") {
            write(nout, format_9974), subnam(1, len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS") {
            write(nout, format_9978), subnam(1, len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen2, p2, "QL") {
        //
        //        xQL:  QL factorization
        //
        if (Mlsamen3, c3, "QLS") {
            write(nout, format_9974), subnam(1, len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS") {
            write(nout, format_9978), subnam(1, len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen2, p2, "RQ") {
        //
        //        xRQ:  RQ factorization
        //
        if (Mlsamen3, c3, "RQS") {
            write(nout, format_9974), subnam(1, len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen5, subnam[(2 - 1) + (6 - 1) * ldsubnam], "LATMS") {
            write(nout, format_9978), subnam(1, len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen2, p2, "LU") {
        //
        if (info != infoe && infoe != 0) {
            write(nout, format_9988), subnam(1, len_trim(subnam)), info, infoe, m, n, n5, imat;
        } else {
            write(nout, format_9975), subnam(1, len_trim(subnam)), info, m, n, n5, imat;
        }
        //
    } else if (Mlsamen2, p2, "CH") {
        //
        if (info != infoe && infoe != 0) {
            write(nout, "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> N =',"
                        "i5,', NB =',i4,', type ',i2)"),
                subnam(1, len_trim(subnam)), info, infoe, m, n5, imat;
        } else {
            write(nout, format_9971), subnam(1, len_trim(subnam)), info, m, n5, imat;
        }
        //
    } else {
        //
        //        PrINTEGER a generic message if the path is unknown.
        //
        write(nout, "(' *** Error code from ',a,' =',i5)"), subnam(1, len_trim(subnam)), info;
    }
    //
    //     Description of error message (alphabetical, left to right)
    //
    //     SUBNAM, INFO, FACT, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, FACT, TRANS, N, KL, KU, NRHS, IMAT
    //
    //     SUBNAM, INFO, FACT, TRANS, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, FACT, UPLO, N, KD, NRHS, IMAT
    //
    //     SUBNAM, INFO, FACT, UPLO, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, FACT, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, FACT, TRANS, N, KL, KU, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, FACT, TRANS, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, FACT, UPLO, N, KD, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, FACT, UPLO, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, M, N, KL, KU, NB, IMAT
    //
    //     SUBNAM, INFO, INFOE, M, N, NB, IMAT
    //
    //     SUBNAM, INFO, INFOE, N, IMAT
    //
    //     SUBNAM, INFO, INFOE, N, KL, KU, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, N, NB, IMAT
    //
    //     SUBNAM, INFO, INFOE, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, UPLO, N, IMAT
    //
    //     SUBNAM, INFO, INFOE, UPLO, N, KD, NB, IMAT
    //
    //     SUBNAM, INFO, INFOE, UPLO, N, KD, NRHS, IMAT
    //
    //     SUBNAM, INFO, INFOE, UPLO, N, NB, IMAT
    //
    //     SUBNAM, INFO, INFOE, UPLO, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, M, N, IMAT
    //
    //     SUBNAM, INFO, M, N, KL, KU, IMAT
    //
    //     SUBNAM, INFO, M, N, KL, KU, NB, IMAT
    //
    //     SUBNAM, INFO, M, N, NB, IMAT
    //
    //     SUBNAM, INFO, M, N, NRHS, NB, IMAT
    //
    //     SUBNAM, INFO, N, IMAT
    //
    //     SUBNAM, INFO, N, KL, KU, NRHS, IMAT
    //
    //     SUBNAM, INFO, N, NB, IMAT
    //
    //     SUBNAM, INFO, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, NORM, N, IMAT
    //
    //     SUBNAM, INFO, NORM, N, KL, KU, IMAT
    //
    //     SUBNAM, INFO, NORM, UPLO, DIAG, N, IMAT
    //
    //     SUBNAM, INFO, NORM, UPLO, DIAG, N, KD, IMAT
    //
    //     SUBNAM, INFO, TRANS, M, N, NRHS, NB, IMAT
    //
    //     SUBNAM, INFO, TRANS, N, KL, KU, NRHS, IMAT
    //
    //     SUBNAM, INFO, TRANS, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, UPLO, DIAG, N, IMAT
    //
    //     SUBNAM, INFO, UPLO, DIAG, N, NB, IMAT
    //
    //     SUBNAM, INFO, UPLO, N, IMAT
    //
    //     SUBNAM, INFO, UPLO, N, KD, IMAT
    //
    //     SUBNAM, INFO, UPLO, N, KD, NB, IMAT
    //
    //     SUBNAM, INFO, UPLO, N, KD, NRHS, IMAT
    //
    //     SUBNAM, INFO, UPLO, N, NB, IMAT
    //
    //     SUBNAM, INFO, UPLO, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, UPLO, TRANS, DIAG, N, KD, NRHS, IMAT
    //
    //     SUBNAM, INFO, UPLO, TRANS, DIAG, N, NRHS, IMAT
    //
    //     SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, IMAT
    //
    //     SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, KD, IMAT
    //
    //     Unknown type
    //
    //     What we do next
    //
    //     End of Alaerh
    //
}
