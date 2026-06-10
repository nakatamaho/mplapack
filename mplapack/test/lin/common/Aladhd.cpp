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

// Derived from LAPACK routine ALADHD.
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

void Aladhd(INTEGER const iounit, fem::str_cref path) {
    common cmn;
    common_write write(cmn);
    //
    // First line of header
    //
    static const char *format_9999 = "(/,1x,a3,' drivers:  General dense matrices')";
    static const char *format_9998 = "(/,1x,a3,' drivers:  General band matrices')";
    static const char *format_9997 = "(/,1x,a3,' drivers:  General tridiagonal')";
    static const char *format_9996 = "(/,1x,a3,' drivers:  ',a9,' positive definite matrices')";
    static const char *format_9995 = "(/,1x,a3,' drivers:  ',a9,' positive definite packed matrices')";
    static const char *format_9994 = "(/,1x,a3,' drivers:  ',a9,' positive definite band matrices')";
    static const char *format_9993 = "(/,1x,a3,' drivers:  ',a9,' positive definite tridiagonal')";
    static const char *format_9971 = "(/,1x,a3,' drivers:  ',a9,' indefinite matrices',', \"Aasen\" Algorithm')";
    static const char *format_9992 = "(/,1x,a3,' drivers:  ',a9,' indefinite matrices',"
                                     "', \"rook\" (bounded Bunch-Kaufman) pivoting')";
    static const char *format_9991 = "(/,1x,a3,' drivers:  ',a9,' indefinite packed matrices',"
                                     "', partial (Bunch-Kaufman) pivoting')";
    static const char *format_9990 = "(/,1x,a3,':  No header available')";
    //
    // GE matrix types
    //
    static const char *format_9989 = "(4x,'1. Diagonal',24x,'7. Last n/2 columns zero',/,4x,"
                                     "'2. Upper triangular',16x,'8. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. Lower triangular',16x,'9. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Random, CNDNUM = 2',13x,'10. Scaled near underflow',/,4x,"
                                     "'5. First column zero',14x,'11. Scaled near overflow',/,4x,"
                                     "'6. Last column zero')";
    //
    // GB matrix types
    //
    static const char *format_9988 = "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. First column zero',15x,'6. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. Last column zero',16x,'7. Scaled near underflow',/,4x,"
                                     "'4. Last n/2 columns zero',11x,'8. Scaled near overflow')";
    //
    // GT matrix types
    //
    static const char *format_9987 = "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                                     "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. First column zero',/,4x,"
                                     "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last column zero',/,4x,"
                                     "'4. Random, CNDNUM = 0.1/EPS',7x,'10. Last n/2 columns zero',/,4x,"
                                     "'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                                     "'6. Scaled near overflow',11x,'12. Scaled near overflow')";
    //
    // PT matrix types
    //
    static const char *format_9986 = "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                                     "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. First row and column zero',/,4x,"
                                     "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last row and column zero',/,"
                                     "4x,'4. Random, CNDNUM = 0.1/EPS',7x,'10. Middle row and column zero',/,"
                                     "4x,'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                                     "'6. Scaled near overflow',11x,'12. Scaled near overflow')";
    //
    // PO, PP matrix types
    //
    static const char *format_9985 = "(4x,'1. Diagonal',24x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = 0.1/EPS',/,3x,"
                                     "'*3. First row and column zero',7x,'8. Scaled near underflow',/,3x,"
                                     "'*4. Last row and column zero',8x,'9. Scaled near overflow',/,3x,"
                                     "'*5. Middle row and column zero',/,3x,'(* - tests error exits from ',a3,"
                                     "'TRF, no test ratios are computed)')";
    //
    // PB matrix types
    //
    static const char *format_9984 = "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,3x,"
                                     "'*2. First row and column zero',7x,'6. Random, CNDNUM = 0.1/EPS',/,3x,"
                                     "'*3. Last row and column zero',8x,'7. Scaled near underflow',/,3x,"
                                     "'*4. Middle row and column zero',6x,'8. Scaled near overflow',/,3x,"
                                     "'(* - tests error exits from ',a3,'TRF, no test ratios are computed)')";
    //
    // SSY, SSP, CHE, CHP matrix types
    //
    static const char *format_9983 = "(4x,'1. Diagonal',24x,'6. Last n/2 rows and columns zero',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. First row and column zero',7x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Last row and column zero',8x,'9. Scaled near underflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'10. Scaled near overflow')";
    //
    // CSY, CSP matrix types
    //
    static const char *format_9982 = "(4x,'1. Diagonal',24x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. First row and column zero',7x,'9. Scaled near underflow',/,4x,"
                                     "'4. Last row and column zero',7x,'10. Scaled near overflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'11. Block diagonal matrix',/,4x,"
                                     "'6. Last n/2 rows and columns zero')";
    //
    // Test ratios
    //
    static const char *format_9981 = "(3x,i2,': norm( L * U - A )  / ( N * norm(A) * EPS )')";
    static const char *format_9980 = "(3x,i2,': norm( B - A * X )  / ','( norm(A) * norm(X) * EPS )')";
    static const char *format_9979 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * CNDNUM * EPS )')";
    static const char *format_9978 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * (error bound) )')";
    static const char *format_9977 = "(3x,i2,': (backward error)   / EPS')";
    static const char *format_9976 = "(3x,i2,': RCOND * CNDNUM - 1.0')";
    static const char *format_9975 = "(3x,i2,': norm( U'' * U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L * L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9974 = "(3x,i2,': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9973 = "(3x,i2,': norm( U''*D*U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9972 = "(3x,i2,': abs( WORK(1) - RPVGRW ) /',' ( max( WORK(1), RPVGRW ) * EPS )')";
    //
    if (iounit <= 0) {
        return;
    }
    fem::str<1> c1 = path(1, 1);
    fem::str<1> c3 = path(3, 3);
    fem::str<2> p2 = path(2, 3);
    bool sord = Mlsame(c1.elems, "S") || Mlsame(c1.elems, "D");
    bool corz = Mlsame(c1.elems, "C") || Mlsame(c1.elems, "Z");
    if (!(sord || corz)) {
        return;
    }
    //
    fem::str<9> sym;
    if (Mlsamen(2, p2.elems, "GE")) {
        //
        // GE: General dense
        //
        write(iounit, format_9999), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9989);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9981), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9978), 4;
        write(iounit, format_9977), 5;
        write(iounit, format_9976), 6;
        write(iounit, format_9972), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "GB")) {
        //
        // GB: General band
        //
        write(iounit, format_9998), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9988);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9981), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9978), 4;
        write(iounit, format_9977), 5;
        write(iounit, format_9976), 6;
        write(iounit, format_9972), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "GT")) {
        //
        // GT: General tridiagonal
        //
        write(iounit, format_9997), path;
        write(iounit, format_9987);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9981), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9978), 4;
        write(iounit, format_9977), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "PO") || Mlsamen(2, p2.elems, "PP") || Mlsamen(2, p2.elems, "PS")) {
        //
        // PO: Positive definite full
        // PS: Positive definite full
        // PP: Positive definite packed
        //
        if (sord) {
            sym = "Symmetric";
        } else {
            sym = "Hermitian";
        }
        if (Mlsame(c3.elems, "O")) {
            write(iounit, format_9996), path, sym;
        } else {
            write(iounit, format_9995), path, sym;
        }
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9985), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9975), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9978), 4;
        write(iounit, format_9977), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "PB")) {
        //
        // PB: Positive definite band
        //
        if (sord) {
            write(iounit, format_9994), path, "Symmetric";
        } else {
            write(iounit, format_9994), path, "Hermitian";
        }
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9984), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9975), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9978), 4;
        write(iounit, format_9977), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "PT")) {
        //
        // PT: Positive definite tridiagonal
        //
        if (sord) {
            write(iounit, format_9993), path, "Symmetric";
        } else {
            write(iounit, format_9993), path, "Hermitian";
        }
        write(iounit, format_9986);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9973), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9978), 4;
        write(iounit, format_9977), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "SY") || Mlsamen(2, p2.elems, "SP")) {
        //
        // SY: Symmetric indefinite full
        // with partial (Bunch-Kaufman) pivoting algorithm
        // SP: Symmetric indefinite packed
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3.elems, "Y")) {
            write(iounit, format_9992), path, "Symmetric";
        } else {
            write(iounit, format_9991), path, "Symmetric";
        }
        write(iounit, "(' Matrix types:')");
        if (sord) {
            write(iounit, format_9983);
        } else {
            write(iounit, format_9982);
        }
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9974), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9977), 4;
        write(iounit, format_9978), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "SR") || Mlsamen(2, p2.elems, "SK")) {
        //
        // SR: Symmetric indefinite full,
        // with rook (bounded Bunch-Kaufman) pivoting algorithm
        //
        // SK: Symmetric indefinite full,
        // with rook (bounded Bunch-Kaufman) pivoting algorithm,
        // ( new storage format for factors:
        // L and diagonal of D is stored in A,
        // subdiagonal of D is stored in E )
        //
        write(iounit, format_9992), path, "Symmetric";
        //
        write(iounit, "(' Matrix types:')");
        if (sord) {
            write(iounit, format_9983);
        } else {
            write(iounit, format_9982);
        }
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9974), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "HA")) {
        //
        // HA: Hermitian
        // Aasen algorithm
        write(iounit, format_9971), path, "Hermitian";
        //
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9983);
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9974), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9977), 4;
        write(iounit, format_9978), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "HE") || Mlsamen(2, p2.elems, "HP")) {
        //
        // HE: Hermitian indefinite full
        // with partial (Bunch-Kaufman) pivoting algorithm
        // HP: Hermitian indefinite packed
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3.elems, "E")) {
            write(iounit, format_9992), path, "Hermitian";
        } else {
            write(iounit, format_9991), path, "Hermitian";
        }
        //
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9983);
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9974), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, format_9977), 4;
        write(iounit, format_9978), 5;
        write(iounit, format_9976), 6;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "HR") || Mlsamen(2, p2.elems, "HK")) {
        //
        // HR: Hermitian indefinite full,
        // with rook (bounded Bunch-Kaufman) pivoting algorithm
        //
        // HK: Hermitian indefinite full,
        // with rook (bounded Bunch-Kaufman) pivoting algorithm,
        // ( new storage format for factors:
        // L and diagonal of D is stored in A,
        // subdiagonal of D is stored in E )
        //
        write(iounit, format_9992), path, "Hermitian";
        //
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9983);
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9974), 1;
        write(iounit, format_9980), 2;
        write(iounit, format_9979), 3;
        write(iounit, "(' Messages:')");
        //
    } else {
        //
        // Print error message if no header is available.
        //
        write(iounit, format_9990), path;
    }
    //
    // End of Aladhd
    //
}
