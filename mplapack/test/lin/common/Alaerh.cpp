/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ALAERH.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Alaerh(fem::str_cref path, fem::str_cref subnam, INTEGER const info, INTEGER const infoe, fem::str_cref opts, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const n5, INTEGER const imat, INTEGER const nfail, INTEGER &nerrs, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    // Description of error message (alphabetical, left to right)
    //
    // SUBNAM, INFO, FACT, N, NRHS, IMAT
    //
    static const char *format_9999 = "(' *** Error code from ',a,'=',i5,', FACT=''',a1,''', N=',i5,', NRHS=',"
                                     "i4,', type ',i2)";
    //
    // SUBNAM, INFO, FACT, TRANS, N, KL, KU, NRHS, IMAT
    //
    static const char *format_9998 = "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,''', TRANS=''',a1,"
                                     "''', N=',i5,', KL=',i5,', KU=',i5,', NRHS=',i4,', type ',i1)";
    //
    // SUBNAM, INFO, FACT, TRANS, N, NRHS, IMAT
    //
    static const char *format_9997 = "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,''', TRANS=''',a1,"
                                     "''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, FACT, UPLO, N, KD, NRHS, IMAT
    //
    static const char *format_9996 = "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,''', UPLO=''',a1,"
                                     "''', N=',i5,', KD=',i5,', NRHS=',i4,', type ',i2)";
    //
    // SUBNAM, INFO, FACT, UPLO, N, NRHS, IMAT
    //
    static const char *format_9995 = "(' *** Error code from ',a,' =',i5,/,' ==> FACT=''',a1,''', UPLO=''',a1,"
                                     "''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, FACT, N, NRHS, IMAT
    //
    static const char *format_9994 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, FACT, TRANS, N, KL, KU, NRHS, IMAT
    //
    static const char *format_9993 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', TRANS=''',a1,''', N=',i5,', KL=',i5,', KU=',i5,', NRHS=',i4,"
                                     "', type ',i1)";
    //
    // SUBNAM, INFO, INFOE, FACT, TRANS, N, NRHS, IMAT
    //
    static const char *format_9992 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', TRANS=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, FACT, UPLO, N, KD, NRHS, IMAT
    //
    static const char *format_9991 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', UPLO=''',a1,''', N=',i5,', KD=',i5,', NRHS=',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, FACT, UPLO, N, NRHS, IMAT
    //
    static const char *format_9990 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> FACT=''',"
                                     "a1,''', UPLO=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, M, N, KL, KU, NB, IMAT
    //
    static const char *format_9989 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> M = ',i5,"
                                     "', N =',i5,', KL =',i5,', KU =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, M, N, NB, IMAT
    //
    static const char *format_9988 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> M =',i5,"
                                     "', N =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, N, IMAT
    //
    static const char *format_9987 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,' for N=',i5,"
                                     "', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, N, KL, KU, NRHS, IMAT
    //
    static const char *format_9986 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> N =',i5,"
                                     "', KL =',i5,', KU =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, N, NB, IMAT
    //
    static const char *format_9985 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> N =',i5,"
                                     "', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, N, NRHS, IMAT
    //
    static const char *format_9984 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> N =',i5,"
                                     "', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, UPLO, N, IMAT
    //
    static const char *format_9983 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                                     "' ==> UPLO = ''',a1,''', N =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, UPLO, N, KD, NB, IMAT
    //
    static const char *format_9982 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                                     "' ==> UPLO = ''',a1,''', N =',i5,', KD =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, UPLO, N, KD, NRHS, IMAT
    //
    static const char *format_9981 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,' ==> UPLO=''',"
                                     "a1,''', N =',i5,', KD =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, UPLO, N, NB, IMAT
    //
    static const char *format_9980 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                                     "' ==> UPLO = ''',a1,''', N =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, INFOE, UPLO, N, NRHS, IMAT
    //
    static const char *format_9979 = "(' *** ',a,' returned with INFO =',i5,' instead of ',i2,/,"
                                     "' ==> UPLO = ''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, M, N, IMAT
    //
    static const char *format_9978 = "(' *** Error code from ',a,' =',i5,' for M =',i5,', N =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, M, N, KL, KU, IMAT
    //
    static const char *format_9977 = "(' *** Error code from ',a,' =',i5,/,' ==> M = ',i5,', N =',i5,', KL =',"
                                     "i5,', KU =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, M, N, KL, KU, NB, IMAT
    //
    static const char *format_9976 = "(' *** Error code from ',a,' =',i5,/,' ==> M = ',i5,', N =',i5,', KL =',"
                                     "i5,', KU =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, M, N, NB, IMAT
    //
    static const char *format_9975 = "(' *** Error code from ',a,'=',i5,' for M=',i5,', N=',i5,', NB=',i4,"
                                     "', type ',i2)";
    //
    // SUBNAM, INFO, M, N, NRHS, NB, IMAT
    //
    static const char *format_9974 = "(' *** Error code from ',a,'=',i5,/,' ==> M =',i5,', N =',i5,', NRHS =',"
                                     "i4,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, N, IMAT
    //
    static const char *format_9973 = "(' *** Error code from ',a,' =',i5,' for N =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, N, KL, KU, NRHS, IMAT
    //
    static const char *format_9972 = "(' *** Error code from ',a,' =',i5,/,' ==> N =',i5,', KL =',i5,', KU =',"
                                     "i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, N, NB, IMAT
    //
    static const char *format_9971 = "(' *** Error code from ',a,'=',i5,' for N=',i5,', NB=',i4,', type ',i2)";
    //
    // SUBNAM, INFO, N, NRHS, IMAT
    //
    static const char *format_9970 = "(' *** Error code from ',a,' =',i5,' for N =',i5,', NRHS =',i4,', type ',"
                                     "i2)";
    //
    // SUBNAM, INFO, NORM, N, IMAT
    //
    static const char *format_9969 = "(' *** Error code from ',a,' =',i5,' for NORM = ''',a1,''', N =',i5,"
                                     "', type ',i2)";
    //
    // SUBNAM, INFO, NORM, N, KL, KU, IMAT
    //
    static const char *format_9968 = "(' *** Error code from ',a,' =',i5,/,' ==> NORM =''',a1,''', N =',i5,"
                                     "', KL =',i5,', KU =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, NORM, UPLO, DIAG, N, IMAT
    //
    static const char *format_9967 = "(' *** Error code from ',a,' =',i5,/,' ==> NORM=''',a1,''', UPLO =''',a1,"
                                     "''', DIAG=''',a1,''', N =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, NORM, UPLO, DIAG, N, KD, IMAT
    //
    static const char *format_9966 = "(' *** Error code from ',a,' =',i5,/,' ==> NORM=''',a1,''', UPLO =''',a1,"
                                     "''', DIAG=''',a1,''', N=',i5,', KD=',i5,', type ',i2)";
    //
    // SUBNAM, INFO, TRANS, M, N, NRHS, NB, IMAT
    //
    static const char *format_9965 = "(' *** Error code from ',a,' =',i5,/,' ==> TRANS = ''',a1,''', M =',i5,"
                                     "', N =',i5,', NRHS =',i4,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, TRANS, N, KL, KU, NRHS, IMAT
    //
    static const char *format_9964 = "(' *** Error code from ',a,'=',i5,/,' ==> TRANS=''',a1,''', N =',i5,"
                                     "', KL =',i5,', KU =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, TRANS, N, NRHS, IMAT
    //
    static const char *format_9963 = "(' *** Error code from ',a,' =',i5,/,' ==> TRANS = ''',a1,''', N =',i5,"
                                     "', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, DIAG, N, IMAT
    //
    static const char *format_9962 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', DIAG =''',a1,"
                                     "''', N =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, DIAG, N, NB, IMAT
    //
    static const char *format_9961 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', DIAG =''',a1,"
                                     "''', N =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, N, IMAT
    //
    static const char *format_9960 = "(' *** Error code from ',a,' =',i5,' for UPLO = ''',a1,''', N =',i5,"
                                     "', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, N, KD, IMAT
    //
    static const char *format_9959 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', KD =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, N, KD, NB, IMAT
    //
    static const char *format_9958 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', KD =',i5,', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, N, KD, NRHS, IMAT
    //
    static const char *format_9957 = "(' *** Error code from ',a,'=',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', KD =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, N, NB, IMAT
    //
    static const char *format_9956 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', NB =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, N, NRHS, IMAT
    //
    static const char *format_9955 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO = ''',a1,''', N =',i5,"
                                     "', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, TRANS, DIAG, N, KD, NRHS, IMAT
    //
    static const char *format_9954 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', TRANS=''',a1,"
                                     "''', DIAG=''',a1,''', N=',i5,', KD=',i5,', NRHS=',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, TRANS, DIAG, N, NRHS, IMAT
    //
    static const char *format_9953 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', TRANS=''',a1,"
                                     "''', DIAG=''',a1,''', N =',i5,', NRHS =',i4,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, IMAT
    //
    static const char *format_9952 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', TRANS=''',a1,"
                                     "''', DIAG=''',a1,''', NORMIN=''',a1,''', N =',i5,', type ',i2)";
    //
    // SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, KD, IMAT
    //
    static const char *format_9951 = "(' *** Error code from ',a,' =',i5,/,' ==> UPLO=''',a1,''', TRANS=''',a1,"
                                     "''', DIAG=''',a1,''', NORMIN=''',a1,''', N=',i5,', KD=',i5,', type ',i2)";
    //
    // Unknown type
    //
    static const char *format_9950 = "(' *** Error code from ',a,' =',i5)";
    //
    // What we do next
    //
    static const char *format_9949 = "(' ==> Doing only the condition estimate for this case')";
    //
    // SUBNAM, INFO, M, N, NB, IMAT
    //
    static const char *format_9930 = "(' *** Error code from ',a,'=',i5,/,' ==> M =',i5,', N =',i5,', NX =',i5,"
                                     "', NB =',i4,', type ',i2)";
    //
    if (info == 0) {
        return;
    }
    fem::str<2> p2 = path(2, 3);
    fem::str<3> c3 = subnam(4, 6);
    //
    // Print the header if this is the first error message.
    //
    if (nfail == 0 && nerrs == 0) {
        if (Mlsamen(3, c3.elems, "SV ") || Mlsamen(3, c3.elems, "SVX")) {
            Aladhd(nout, path);
        } else {
            Alahd(nout, path);
        }
    }
    nerrs++;
    //
    // Print the message detailing the error and form of recovery,
    // if any.
    //
    fem::str<1> uplo;
    if (Mlsamen(2, p2.elems, "GE")) {
        //
        // xGE:  General matrices
        //
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9988), subnam(1, fem::len_trim(subnam)), info, infoe, m, n, n5, imat;
            } else {
                write(nout, format_9975), subnam(1, fem::len_trim(subnam)), info, m, n, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9984), subnam(1, fem::len_trim(subnam)), info, infoe, n, n5, imat;
            } else {
                write(nout, format_9970), subnam(1, fem::len_trim(subnam)), info, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9992), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9997), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "TRI")) {
            //
            write(nout, format_9971), subnam(1, fem::len_trim(subnam)), info, n, n5, imat;
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            //
            write(nout, format_9978), subnam(1, fem::len_trim(subnam)), info, m, n, imat;
            //
        } else if (Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9969), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, imat;
            //
        } else if (Mlsamen(3, c3.elems, "LS ")) {
            //
            write(nout, format_9965), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, n, kl, n5, imat;
            //
        } else if (Mlsamen(3, c3.elems, "LSX") || Mlsamen(3, c3.elems, "LSS")) {
            //
            write(nout, format_9974), subnam(1, fem::len_trim(subnam)), info, m, n, kl, n5, imat;
            //
        } else {
            //
            write(nout, format_9963), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "GB")) {
        //
        // xGB:  General band matrices
        //
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9989), subnam(1, fem::len_trim(subnam)), info, infoe, m, n, kl, ku, n5, imat;
            } else {
                write(nout, format_9976), subnam(1, fem::len_trim(subnam)), info, m, n, kl, ku, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9986), subnam(1, fem::len_trim(subnam)), info, infoe, n, kl, ku, n5, imat;
            } else {
                write(nout, format_9972), subnam(1, fem::len_trim(subnam)), info, n, kl, ku, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9993), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, kl, ku, n5, imat;
            } else {
                write(nout, format_9998), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, kl, ku, n5, imat;
            }
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            //
            write(nout, format_9977), subnam(1, fem::len_trim(subnam)), info, m, n, kl, ku, imat;
            //
        } else if (Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9968), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, kl, ku, imat;
            //
        } else {
            //
            write(nout, format_9964), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, kl, ku, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "GT")) {
        //
        // xGT:  General tridiagonal matrices
        //
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9987), subnam(1, fem::len_trim(subnam)), info, infoe, n, imat;
            } else {
                write(nout, format_9973), subnam(1, fem::len_trim(subnam)), info, n, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9984), subnam(1, fem::len_trim(subnam)), info, infoe, n, n5, imat;
            } else {
                write(nout, format_9970), subnam(1, fem::len_trim(subnam)), info, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9992), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9997), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9969), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, imat;
            //
        } else {
            //
            write(nout, format_9963), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "PO")) {
        //
        // xPO:  Symmetric or Hermitian positive definite matrices
        //
        uplo = opts(1, 1);
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9980), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, m, n5, imat;
            } else {
                write(nout, format_9956), subnam(1, fem::len_trim(subnam)), info, uplo, m, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam(1, fem::len_trim(subnam)), info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "TRI")) {
            //
            write(nout, format_9956), subnam(1, fem::len_trim(subnam)), info, uplo, m, n5, imat;
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS") || Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9960), subnam(1, fem::len_trim(subnam)), info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam(1, fem::len_trim(subnam)), info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "PS")) {
        //
        // xPS:  Symmetric or Hermitian positive semi-definite matrices
        //
        uplo = opts(1, 1);
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9980), subnam, info, infoe, uplo, m, n5, imat;
            } else {
                write(nout, format_9956), subnam, info, uplo, m, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam, info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam, info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam, info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam, info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "TRI")) {
            //
            write(nout, format_9956), subnam, info, uplo, m, n5, imat;
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMT") || Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9960), subnam, info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam, info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "SY") || Mlsamen(2, p2.elems, "SR") || Mlsamen(2, p2.elems, "SK") || Mlsamen(2, p2.elems, "HE") || Mlsamen(2, p2.elems, "HR") || Mlsamen(2, p2.elems, "HK") || Mlsamen(2, p2.elems, "HA")) {
        //
        // xSY: symmetric indefinite matrices
        // with partial (Bunch-Kaufman) pivoting;
        // xSR: symmetric indefinite matrices
        // with rook (bounded Bunch-Kaufman) pivoting;
        // xSK: symmetric indefinite matrices
        // with rook (bounded Bunch-Kaufman) pivoting,
        // new storage format;
        // xHE: Hermitian indefinite matrices
        // with partial (Bunch-Kaufman) pivoting.
        // xHR: Hermitian indefinite matrices
        // with rook (bounded Bunch-Kaufman) pivoting;
        // xHK: Hermitian indefinite matrices
        // with rook (bounded Bunch-Kaufman) pivoting,
        // new storage format;
        // xHA: Hermitian matrices
        // Aasen Algorithm
        //
        uplo = opts(1, 1);
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9980), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, m, n5, imat;
            } else {
                write(nout, format_9956), subnam(1, fem::len_trim(subnam)), info, uplo, m, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(2, c3.elems, "SV")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam(1, fem::len_trim(subnam)), info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS") || Mlsamen(3, c3.elems, "TRI") || Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9960), subnam(1, fem::len_trim(subnam)), info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam(1, fem::len_trim(subnam)), info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "PP") || Mlsamen(2, p2.elems, "SP") || Mlsamen(2, p2.elems, "HP")) {
        //
        // xPP, xHP, or xSP:  Symmetric or Hermitian packed matrices
        //
        uplo = opts(1, 1);
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9983), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, m, imat;
            } else {
                write(nout, format_9960), subnam(1, fem::len_trim(subnam)), info, uplo, m, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9979), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, n, n5, imat;
            } else {
                write(nout, format_9955), subnam(1, fem::len_trim(subnam)), info, uplo, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9990), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, n5, imat;
            } else {
                write(nout, format_9995), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, n5, imat;
            }
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS") || Mlsamen(3, c3.elems, "TRI") || Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9960), subnam(1, fem::len_trim(subnam)), info, uplo, m, imat;
            //
        } else {
            //
            write(nout, format_9955), subnam(1, fem::len_trim(subnam)), info, uplo, m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "PB")) {
        //
        // xPB:  Symmetric (Hermitian) positive definite band matrix
        //
        uplo = opts(1, 1);
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9982), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, m, kl, n5, imat;
            } else {
                write(nout, format_9958), subnam(1, fem::len_trim(subnam)), info, uplo, m, kl, n5, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9981), subnam(1, fem::len_trim(subnam)), info, infoe, uplo, n, kl, n5, imat;
            } else {
                write(nout, format_9957), subnam(1, fem::len_trim(subnam)), info, uplo, n, kl, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9991), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), opts(2, 2), n, kl, n5, imat;
            } else {
                write(nout, format_9996), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), n, kl, n5, imat;
            }
            //
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS") || Mlsamen(3, c3.elems, "CON")) {
            //
            write(nout, format_9959), subnam(1, fem::len_trim(subnam)), info, uplo, m, kl, imat;
            //
        } else {
            //
            write(nout, format_9957), subnam(1, fem::len_trim(subnam)), info, uplo, m, kl, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "PT")) {
        //
        // xPT:  Positive definite tridiagonal matrices
        //
        if (Mlsamen(3, c3.elems, "TRF")) {
            if (info != infoe && infoe != 0) {
                write(nout, format_9987), subnam(1, fem::len_trim(subnam)), info, infoe, n, imat;
            } else {
                write(nout, format_9973), subnam(1, fem::len_trim(subnam)), info, n, imat;
            }
            if (info != 0) {
                write(nout, format_9949);
            }
            //
        } else if (Mlsamen(3, c3.elems, "SV ")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9984), subnam(1, fem::len_trim(subnam)), info, infoe, n, n5, imat;
            } else {
                write(nout, format_9970), subnam(1, fem::len_trim(subnam)), info, n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "SVX")) {
            //
            if (info != infoe && infoe != 0) {
                write(nout, format_9994), subnam(1, fem::len_trim(subnam)), info, infoe, opts(1, 1), n, n5, imat;
            } else {
                write(nout, format_9999), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), n, n5, imat;
            }
            //
        } else if (Mlsamen(3, c3.elems, "CON")) {
            //
            if (Mlsame(subnam(1, 1).elems(), "S") || Mlsame(subnam(1, 1).elems(), "D")) {
                write(nout, format_9973), subnam(1, fem::len_trim(subnam)), info, m, imat;
            } else {
                write(nout, format_9969), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, imat;
            }
            //
        } else {
            //
            write(nout, format_9963), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "TR")) {
        //
        // xTR:  Triangular matrix
        //
        if (Mlsamen(3, c3.elems, "TRI")) {
            write(nout, format_9961), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), m, n5, imat;
        } else if (Mlsamen(3, c3.elems, "CON")) {
            write(nout, format_9967), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATRS")) {
            write(nout, format_9952), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), opts(4, 4), m, imat;
        } else {
            write(nout, format_9953), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "TP")) {
        //
        // xTP:  Triangular packed matrix
        //
        if (Mlsamen(3, c3.elems, "TRI")) {
            write(nout, format_9962), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), m, imat;
        } else if (Mlsamen(3, c3.elems, "CON")) {
            write(nout, format_9967), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATPS")) {
            write(nout, format_9952), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), opts(4, 4), m, imat;
        } else {
            write(nout, format_9953), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "TB")) {
        //
        // xTB:  Triangular band matrix
        //
        if (Mlsamen(3, c3.elems, "CON")) {
            write(nout, format_9966), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, kl, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATBS")) {
            write(nout, format_9951), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), opts(4, 4), m, kl, imat;
        } else {
            write(nout, format_9954), subnam(1, fem::len_trim(subnam)), info, opts(1, 1), opts(2, 2), opts(3, 3), m, kl, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "QR")) {
        //
        // xQR:  QR factorization
        //
        if (Mlsamen(3, c3.elems, "QRS")) {
            write(nout, format_9974), subnam(1, fem::len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            write(nout, format_9978), subnam(1, fem::len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "QK")) {
        //
        // xQK:  truncated QR factorization with pivoting
        //
        if (Mlsamen(7, subnam(2, 8).elems(), "GEQP3RK")) {
            write(nout, format_9930), subnam(1, fem::len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            write(nout, format_9978), subnam(1, fem::len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "LQ")) {
        //
        // xLQ:  LQ factorization
        //
        if (Mlsamen(3, c3.elems, "LQS")) {
            write(nout, format_9974), subnam(1, fem::len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            write(nout, format_9978), subnam(1, fem::len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "QL")) {
        //
        // xQL:  QL factorization
        //
        if (Mlsamen(3, c3.elems, "QLS")) {
            write(nout, format_9974), subnam(1, fem::len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            write(nout, format_9978), subnam(1, fem::len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "RQ")) {
        //
        // xRQ:  RQ factorization
        //
        if (Mlsamen(3, c3.elems, "RQS")) {
            write(nout, format_9974), subnam(1, fem::len_trim(subnam)), info, m, n, kl, n5, imat;
        } else if (Mlsamen(5, subnam(2, 6).elems(), "LATMS")) {
            write(nout, format_9978), subnam(1, fem::len_trim(subnam)), info, m, n, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "LU")) {
        //
        if (info != infoe && infoe != 0) {
            write(nout, format_9988), subnam(1, fem::len_trim(subnam)), info, infoe, m, n, n5, imat;
        } else {
            write(nout, format_9975), subnam(1, fem::len_trim(subnam)), info, m, n, n5, imat;
        }
        //
    } else if (Mlsamen(2, p2.elems, "CH")) {
        //
        if (info != infoe && infoe != 0) {
            write(nout, format_9985), subnam(1, fem::len_trim(subnam)), info, infoe, m, n5, imat;
        } else {
            write(nout, format_9971), subnam(1, fem::len_trim(subnam)), info, m, n5, imat;
        }
        //
    } else {
        //
        // Print a generic message if the path is unknown.
        //
        write(nout, format_9950), subnam(1, fem::len_trim(subnam)), info;
    }
    //
    // End of Alaerh
    //
}
