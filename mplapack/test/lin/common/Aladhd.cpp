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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

void Aladhd(common &cmn, INTEGER const iounit, const char *path) {
    common_write write(cmn);
    static const char *format_9972 = "(3x,i2,': abs( WORK(1) - RPVGRW ) /',' ( max( WORK(1), RPVGRW ) * EPS )')";
    static const char *format_9974 = "(3x,i2,': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9975 = "(3x,i2,': norm( U'' * U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L * L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9976 = "(3x,i2,': RCOND * CNDNUM - 1.0')";
    static const char *format_9977 = "(3x,i2,': (backward error)   / EPS')";
    static const char *format_9978 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * (error bound) )')";
    static const char *format_9979 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * CNDNUM * EPS )')";
    static const char *format_9980 = "(3x,i2,': norm( B - A * X )  / ','( norm(A) * norm(X) * EPS )')";
    static const char *format_9981 = "(3x,i2,': norm( L * U - A )  / ( N * norm(A) * EPS )')";
    static const char *format_9982 = "(4x,'1. Diagonal',24x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. First row and column zero',7x,'9. Scaled near underflow',/,4x,"
                                     "'4. Last row and column zero',7x,'10. Scaled near overflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'11. Block diagonal matrix',/,4x,"
                                     "'6. Last n/2 rows and columns zero')";
    static const char *format_9983 = "(4x,'1. Diagonal',24x,'6. Last n/2 rows and columns zero',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. First row and column zero',7x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Last row and column zero',8x,'9. Scaled near underflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'10. Scaled near overflow')";
    static const char *format_9991 = "(/,1x,a3,' drivers:  ',a9,' indefinite packed matrices',"
                                     "', partial (Bunch-Kaufman) pivoting')";
    static const char *format_9992 = "(/,1x,a3,' drivers:  ',a9,' indefinite matrices',"
                                     "', \"rook\" (bounded Bunch-Kaufman) pivoting')";
    static const char *format_9993 = "(/,1x,a3,' drivers:  ',a9,' positive definite tridiagonal')";
    static const char *format_9994 = "(/,1x,a3,' drivers:  ',a9,' positive definite band matrices')";
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
    //     .. Executable Statements ..
    //
    if (iounit <= 0) {
        return;
    }
    char c1[1];
    c1[0] = path[0];
    char c3[1];
    c3[0] = path[2];
    char p2[2];
    p2[0] = path[1];
    p2[1] = path[2];
    bool sord = Mlsame(c1, "S") || Mlsame(c1, "D");
    bool corz = Mlsame(c1, "C") || Mlsame(c1, "Z");
    if (!(sord || corz)) {
        return;
    }
    //
    char sym[9];
    if (Mlsamen(2, p2, "GE")) {
        //
        //        GE: General dense
        //
        write(6, "(/,1x,a3,' drivers:  General dense matrices')"), path;
        write(6, "(' Matrix types:')");
        write(6, "(4x,'1. Diagonal',24x,'7. Last n/2 columns zero',/,4x,"
                 "'2. Upper triangular',16x,'8. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                 "'3. Lower triangular',16x,'9. Random, CNDNUM = 0.1/EPS',/,4x,"
                 "'4. Random, CNDNUM = 2',13x,'10. Scaled near underflow',/,4x,"
                 "'5. First column zero',14x,'11. Scaled near overflow',/,4x,"
                 "'6. Last column zero')");
        write(6, "(' Test ratios:')");
        write(6, format_9981), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9978), 4;
        write(6, format_9977), 5;
        write(6, format_9976), 6;
        write(6, format_9972), 7;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "GB")) {
        //
        //        GB: General band
        //
        write(6, "(/,1x,a3,' drivers:  General band matrices')"), path;
        write(6, "(' Matrix types:')");
        write(6, "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,"
                 "4x,'2. First column zero',15x,'6. Random, CNDNUM = 0.1/EPS',/,4x,"
                 "'3. Last column zero',16x,'7. Scaled near underflow',/,4x,"
                 "'4. Last n/2 columns zero',11x,'8. Scaled near overflow')");
        write(6, "(' Test ratios:')");
        write(6, format_9981), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9978), 4;
        write(6, format_9977), 5;
        write(6, format_9976), 6;
        write(6, format_9972), 7;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "GT")) {
        //
        //        GT: General tridiagonal
        //
        write(6, "(/,1x,a3,' drivers:  General tridiagonal')"), path;
        write(6, "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                 "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                 "'2. Random, CNDNUM = 2',14x,'8. First column zero',/,4x,"
                 "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last column zero',/,4x,"
                 "'4. Random, CNDNUM = 0.1/EPS',7x,'10. Last n/2 columns zero',/,4x,"
                 "'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                 "'6. Scaled near overflow',11x,'12. Scaled near overflow')");
        write(6, "(' Test ratios:')");
        write(6, format_9981), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9978), 4;
        write(6, format_9977), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "PO") || Mlsamen(2, p2, "PP") || Mlsamen(2, p2, "PS")) {
        //
        //        PO: Positive definite full
        //        PS: Positive definite full
        //        PP: Positive definite packed
        //
        if (sord) {
            strncpy(sym, "Symmetric", strlen(sym));
        } else {
            strncpy(sym, "Hermitian", strlen(sym));
        }
        if (Mlsame(c3, "O")) {
            write(6, "(/,1x,a3,' drivers:  ',a9,' positive definite matrices')"), path, sym;
        } else {
            write(6, "(/,1x,a3,' drivers:  ',a9,' positive definite packed matrices')"), path, sym;
        }
        write(6, "(' Matrix types:')");
        write(6, "(4x,'1. Diagonal',24x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                 "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = 0.1/EPS',/,3x,"
                 "'*3. First row and column zero',7x,'8. Scaled near underflow',/,3x,"
                 "'*4. Last row and column zero',8x,'9. Scaled near overflow',/,3x,"
                 "'*5. Middle row and column zero',/,3x,'(* - tests error exits from ',"
                 "a3,'TRF, no test ratios are computed)')"),
            path;
        write(6, "(' Test ratios:')");
        write(6, format_9975), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9978), 4;
        write(6, format_9977), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "PB")) {
        //
        //        PB: Positive definite band
        //
        if (sord) {
            write(6, format_9994), path, "Symmetric";
        } else {
            write(6, format_9994), path, "Hermitian";
        }
        write(6, "(' Matrix types:')");
        write(6, "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,"
                 "3x,'*2. First row and column zero',7x,'6. Random, CNDNUM = 0.1/EPS',/,"
                 "3x,'*3. Last row and column zero',8x,'7. Scaled near underflow',/,3x,"
                 "'*4. Middle row and column zero',6x,'8. Scaled near overflow',/,3x,"
                 "'(* - tests error exits from ',a3,'TRF, no test ratios are computed)')"),
            path;
        write(6, "(' Test ratios:')");
        write(6, format_9975), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9978), 4;
        write(6, format_9977), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "PT")) {
        //
        //        PT: Positive definite tridiagonal
        //
        if (sord) {
            write(6, format_9993), path, "Symmetric";
        } else {
            write(6, format_9993), path, "Hermitian";
        }
        write(6, "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                 "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                 "'2. Random, CNDNUM = 2',14x,'8. First row and column zero',/,4x,"
                 "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last row and column zero',/,"
                 "4x,'4. Random, CNDNUM = 0.1/EPS',7x,'10. Middle row and column zero',/,"
                 "4x,'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                 "'6. Scaled near overflow',11x,'12. Scaled near overflow')");
        write(6, "(' Test ratios:')");
        write(6, "(3x,i2,': norm( U''*D*U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                 "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')"),
            1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9978), 4;
        write(6, format_9977), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "SY") || Mlsamen(2, p2, "SP")) {
        //
        //        SY: Symmetric indefinite full
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //        SP: Symmetric indefinite packed
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3, "Y")) {
            write(6, format_9992), path, "Symmetric";
        } else {
            write(6, format_9991), path, "Symmetric";
        }
        write(6, "(' Matrix types:')");
        if (sord) {
            write(6, format_9983);
        } else {
            write(6, format_9982);
        }
        write(6, "(' Test ratios:')");
        write(6, format_9974), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9977), 4;
        write(6, format_9978), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "SR") || Mlsamen(2, p2, "SK")) {
        //
        //        SR: Symmetric indefinite full,
        //            with rook (bounded Bunch-Kaufman) pivoting algorithm
        //
        //        SK: Symmetric indefinite full,
        //            with rook (bounded Bunch-Kaufman) pivoting algorithm,
        //            ( new storage format for factors:
        //              L and diagonal of D is stored in A,
        //              subdiagonal of D is stored in E )
        //
        write(6, format_9992), path, "Symmetric";
        //
        write(6, "(' Matrix types:')");
        if (sord) {
            write(6, format_9983);
        } else {
            write(6, format_9982);
        }
        //
        write(6, "(' Test ratios:')");
        write(6, format_9974), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HA")) {
        //
        //        HA: Hermitian
        //            Aasen algorithm
        write(6, "(/,1x,a3,' drivers:  ',a9,' indefinite matrices',"
                 "', \"Aasen\" Algorithm')"),
            path, "Hermitian";
        //
        write(6, "(' Matrix types:')");
        write(6, format_9983);
        //
        write(6, "(' Test ratios:')");
        write(6, format_9974), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9977), 4;
        write(6, format_9978), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HE") || Mlsamen(2, p2, "HP")) {
        //
        //        HE: Hermitian indefinite full
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //        HP: Hermitian indefinite packed
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3, "E")) {
            write(6, format_9992), path, "Hermitian";
        } else {
            write(6, format_9991), path, "Hermitian";
        }
        //
        write(6, "(' Matrix types:')");
        write(6, format_9983);
        //
        write(6, "(' Test ratios:')");
        write(6, format_9974), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, format_9977), 4;
        write(6, format_9978), 5;
        write(6, format_9976), 6;
        write(6, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HR") || Mlsamen(2, p2, "HK")) {
        //
        //        HR: Hermitian indefinite full,
        //            with rook (bounded Bunch-Kaufman) pivoting algorithm
        //
        //        HK: Hermitian indefinite full,
        //            with rook (bounded Bunch-Kaufman) pivoting algorithm,
        //            ( new storage format for factors:
        //              L and diagonal of D is stored in A,
        //              subdiagonal of D is stored in E )
        //
        write(6, format_9992), path, "Hermitian";
        //
        write(6, "(' Matrix types:')");
        write(6, format_9983);
        //
        write(6, "(' Test ratios:')");
        write(6, format_9974), 1;
        write(6, format_9980), 2;
        write(6, format_9979), 3;
        write(6, "(' Messages:')");
        //
    } else {
        //
        //        Print error message if no header is available.
        //
        write(6, "(/,1x,a3,':  No header available')"), path;
    }
    //
    //     First line of header
    //
    //     GE matrix types
    //
    //     GB matrix types
    //
    //     GT matrix types
    //
    //     PT matrix types
    //
    //     PO, PP matrix types
    //
    //     PB matrix types
    //
    //     SSY, SSP, CHE, CHP matrix types
    //
    //     CSY, CSP matrix types
    //
    //     Test ratios
    //
    //     End of Aladhd
    //
}
