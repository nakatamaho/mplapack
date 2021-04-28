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
#include <string>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_lin.h>

void Alahd(INTEGER const iounit, const char *path) {
    common cmn;
    common_write write(cmn);
    static const char *format_9892 = "(/,1x,a3,':  ',a9,' indefinite matrices',"
                                     "', \"rook\" (bounded Bunch-Kaufman) pivoting')";
    static const char *format_9926 = "(3x,i2,': Largest 2-Norm of 2-by-2 pivots',/,12x,"
                                     "' - ( ( 1 + ALPHA ) / ( 1 - ALPHA ) ) + THRESH')";
    static const char *format_9927 = "(3x,i2,': ABS( Largest element in L )',/,12x,"
                                     "' - ( 1 / ( 1 - ALPHA ) ) + THRESH')";
    static const char *format_9928 = "(7x,'where ALPHA = ( 1 + SQRT( 17 ) ) / 8')";
    static const char *format_9935 = "(3x,i2,': norm( B - A * X )   / ',"
                                     "'( max(M,N) * norm(A) * norm(X) * EPS )')";
    static const char *format_9938 = "(3x,i2,': norm( I - Q''*Q )      / ( M * EPS )')";
    static const char *format_9940 = "(3x,i2,': norm(svd(A) - svd(R)) / ','( M * norm(svd(R)) * EPS )')";
    static const char *format_9941 = "(3x,i2,': norm( C*Q'' - C*Q'' )/ ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9942 = "(3x,i2,': norm( Q''*C - Q''*C )/ ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9943 = "(3x,i2,': norm( C*Q - C*Q )  / ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9944 = "(3x,i2,': norm( Q*C - Q*C )  / ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9945 = "(3x,i2,': norm( I - Q*Q'' )   / ( N * EPS )')";
    static const char *format_9946 = "(3x,i2,': norm( I - Q''*Q )   / ( M * EPS )')";
    static const char *format_9951 = "(' Test ratio for ',a,':',/,3x,i2,"
                                     "': norm( s*b - A*x )  / ( norm(A) * norm(x) * EPS )')";
    static const char *format_9953 = "(3x,i2,': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9954 = "(3x,i2,': norm( U'' * U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L * L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9955 = "(3x,i2,': RCOND * CNDNUM - 1.0')";
    static const char *format_9956 = "(3x,i2,': (backward error)   / EPS')";
    static const char *format_9957 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * (error bound) )')";
    static const char *format_9958 = "(3x,i2,': norm( X - XACT )   / ',"
                                     "'( norm(XACT) * CNDNUM * EPS ), refined')";
    static const char *format_9959 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * CNDNUM * EPS )')";
    static const char *format_9960 = "(3x,i2,': norm( B - A * X )  / ','( norm(A) * norm(X) * EPS )')";
    static const char *format_9961 = "(3x,i2,': norm( I - A*AINV ) / ','( N * norm(A) * norm(AINV) * EPS )')";
    static const char *format_9962 = "(3x,i2,': norm( L * U - A )  / ( N * norm(A) * EPS )')";
    static const char *format_9970 = "(4x,'1. Diagonal',24x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Upper triangular',16x,'6. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. Lower triangular',16x,'7. Scaled near underflow',/,4x,"
                                     "'4. Random, CNDNUM = 2',14x,'8. Scaled near overflow')";
    static const char *format_9971 = "(4x,'1. Diagonal',24x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. First row and column zero',7x,'9. Scaled near underflow',/,4x,"
                                     "'4. Last row and column zero',7x,'10. Scaled near overflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'11. Block diagonal matrix',/,4x,"
                                     "'6. Last n/2 rows and columns zero')";
    static const char *format_9972 = "(4x,'1. Diagonal',24x,'6. Last n/2 rows and columns zero',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. First row and column zero',7x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Last row and column zero',8x,'9. Scaled near underflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'10. Scaled near overflow')";
    static const char *format_9979 = "(4x,'1. Diagonal',24x,'7. Last n/2 columns zero',/,4x,"
                                     "'2. Upper triangular',16x,'8. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. Lower triangular',16x,'9. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Random, CNDNUM = 2',13x,'10. Scaled near underflow',/,4x,"
                                     "'5. First column zero',14x,'11. Scaled near overflow',/,4x,"
                                     "'6. Last column zero')";
    static const char *format_9987 = "(/,1x,a3,':  ',a2,' factorization of general matrices')";
    static const char *format_9991 = "(/,1x,a3,':  ',a9,' indefinite packed matrices',"
                                     "', partial (Bunch-Kaufman) pivoting')";
    static const char *format_9992 = "(/,1x,a3,':  ',a9,' indefinite matrices',"
                                     "', partial (Bunch-Kaufman) pivoting')";
    static const char *format_9993 = "(/,1x,a3,':  ',a9,' positive definite tridiagonal')";
    static const char *format_9994 = "(/,1x,a3,':  ',a9,' positive definite band matrices')";
    static const char *format_9995 = "(/,1x,a3,':  ',a9,' positive definite packed matrices')";
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
    //     .. Executable Statements ..
    //
	if (iounit <= 0) {
        return;
    }
    char c1[1];
    char c3[3];
    char p2[2];
    c1[0] = path[0];
    c3[0] = path[2];
    c3[1] = path[3];
    c3[2] = path[4];
    p2[0] = path[1];            
    p2[1] = path[2];
    bool sord = Mlsame(c1, "R");
    bool corz = Mlsame(c1, "C");
    if (!(sord || corz)) {
        return;
    }
    //
    char sym[10];
    char eigcnm[5];
    std::string subnam;
      if (Mlsamen(2, p2, "GE")) {
        //
        //        GE: General dense
        //
        write(iounit, "(/,1x,a3,':  General dense matrices')"), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9979);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9962), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9959), 4;
        write(iounit, format_9958), 5;
        write(iounit, format_9957), 6;
        write(iounit, format_9956), 7;
        write(iounit, format_9955), 8;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "GB")) {
        //
        //        GB: General band
        //
        write(iounit, "(/,1x,a3,':  General band matrices')"), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,"
                      "4x,'2. First column zero',15x,'6. Random, CNDNUM = .01/EPS',/,4x,"
                      "'3. Last column zero',16x,'7. Scaled near underflow',/,4x,"
                      "'4. Last n/2 columns zero',11x,'8. Scaled near overflow')");
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9962), 1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "GT")) {
        //
        //        GT: General tridiagonal
        //
        write(iounit, "(/,1x,a3,':  General tridiagonal')"), path;
        write(iounit, "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                      "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                      "'2. Random, CNDNUM = 2',14x,'8. First column zero',/,4x,"
                      "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last column zero',/,4x,"
                      "'4. Random, CNDNUM = 0.1/EPS',7x,'10. Last n/2 columns zero',/,4x,"
                      "'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                      "'6. Scaled near overflow',11x,'12. Scaled near overflow')");
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9962), 1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "PO") || Mlsamen(2, p2, "PP")) {
        //
        //        PO: Positive definite full
        //        PP: Positive definite packed
        //
        if (sord) {
            strncpy(sym, "Symmetric", strlen(sym) - 1);
        } else {
            strncpy(sym, "Hermitian", strlen(sym) - 1);
        }
        if (Mlsame(c3, "O")) {
            write(iounit, "(/,1x,a3,':  ',a9,' positive definite matrices')"), path, sym;
        } else {
            write(iounit, format_9995), path, sym;
        }
        write(iounit, "(' Matrix types:')");
        write(iounit, "(4x,'1. Diagonal',24x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                      "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = 0.1/EPS',/,3x,"
                      "'*3. First row and column zero',7x,'8. Scaled near underflow',/,3x,"
                      "'*4. Last row and column zero',8x,'9. Scaled near overflow',/,3x,"
                      "'*5. Middle row and column zero',/,3x,'(* - tests error exits from ',"
                      "a3,'TRF, no test ratios are computed)')"),
            path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9954), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9959), 4;
        write(iounit, format_9958), 5;
        write(iounit, format_9957), 6;
        write(iounit, format_9956), 7;
        write(iounit, format_9955), 8;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "PS")) {
        //
        //        PS: Positive semi-definite full
        //
        if (sord) {
            strncpy(sym, "Symmetric", strlen(sym) - 1);
        } else {
            strncpy(sym, "Hermitian", strlen(sym) - 1);
        }
        if (Mlsame(c1, "S") || Mlsame(c1, "C")) {
            strncpy(eigcnm, "1E04", strlen(eigcnm) - 1);
        } else {
            strncpy(eigcnm, "1D12", strlen(eigcnm) - 1);
        }
        write(iounit, format_9995), path, sym;
        write(iounit, "(' Matrix types:')");
        write(iounit, "(4x,'1. Diagonal',/,4x,'2. Random, CNDNUM = 2',14x,/,3x,"
                      "'*3. Nonzero eigenvalues of: D(1:RANK-1)=1 and ','D(RANK) = 1.0/',a4,/,"
                      "3x,'*4. Nonzero eigenvalues of: D(1)=1 and ',' D(2:RANK) = 1.0/',a4,/,"
                      "3x,'*5. Nonzero eigenvalues of: D(I) = ',a4,'**(-(I-1)/(RANK-1)) ',"
                      "' I=1:RANK',/,4x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                      "'7. Random, CNDNUM = 0.1/EPS',/,4x,'8. Scaled near underflow',/,4x,"
                      "'9. Scaled near overflow',/,3x,'(* - Semi-definite tests )')"),
            eigcnm, eigcnm, eigcnm;
        write(iounit, "(' Difference:')");
        write(iounit, "(3x,'RANK minus computed rank, returned by ',a,'PSTRF')"), c1;
        write(iounit, "(' Test ratio:')");
        write(iounit, "(3x,'norm( P * U'' * U * P'' - A ) / ( N * norm(A) * EPS )',', or',/,"
                      "3x,'norm( P * L * L'' * P'' - A ) / ( N * norm(A) * EPS )')");
        write(iounit, "(' Messages:')");
    } else if (Mlsamen(2, p2, "PB")) {
        //
        //        PB: Positive definite band
        //
        if (sord) {
            write(iounit, format_9994), path, "Symmetric";
        } else {
            write(iounit, format_9994), path, "Hermitian";
        }
        write(iounit, "(' Matrix types:')");
        write(iounit, "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,"
                      "3x,'*2. First row and column zero',7x,'6. Random, CNDNUM = 0.1/EPS',/,"
                      "3x,'*3. Last row and column zero',8x,'7. Scaled near underflow',/,3x,"
                      "'*4. Middle row and column zero',6x,'8. Scaled near overflow',/,3x,"
                      "'(* - tests error exits from ',a3,'TRF, no test ratios are computed)')"),
            path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9954), 1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "PT")) {
        //
        //        PT: Positive definite tridiagonal
        //
        if (sord) {
            write(iounit, format_9993), path, "Symmetric";
        } else {
            write(iounit, format_9993), path, "Hermitian";
        }
        write(iounit, "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                      "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                      "'2. Random, CNDNUM = 2',14x,'8. First row and column zero',/,4x,"
                      "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last row and column zero',/,"
                      "4x,'4. Random, CNDNUM = 0.1/EPS',7x,'10. Middle row and column zero',/,"
                      "4x,'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                      "'6. Scaled near overflow',11x,'12. Scaled near overflow')");
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( U''*D*U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                      "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')"),
            1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "SY")) {
        //
        //        SY: Symmetric indefinite full,
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3, "Y")) {
            write(iounit, format_9992), path, "Symmetric";
        } else {
            write(iounit, format_9991), path, "Symmetric";
        }
        write(iounit, "(' Matrix types:')");
        if (sord) {
            write(iounit, format_9972);
        } else {
            write(iounit, format_9971);
        }
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9960), 4;
        write(iounit, format_9959), 5;
        write(iounit, format_9958), 6;
        write(iounit, format_9956), 7;
        write(iounit, format_9957), 8;
        write(iounit, format_9955), 9;
        write(iounit, "(' Messages:')");
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
        write(iounit, format_9892), path, "Symmetric";
        //
        write(iounit, "(' Matrix types:')");
        if (sord) {
            write(iounit, format_9972);
        } else {
            write(iounit, format_9971);
        }
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9927), 3;
        write(iounit, format_9928);
        write(iounit, format_9926), 4;
        write(iounit, format_9928);
        write(iounit, format_9960), 5;
        write(iounit, format_9959), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "SP")) {
        //
        //        SP: Symmetric indefinite packed,
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3, "Y")) {
            write(iounit, format_9992), path, "Symmetric";
        } else {
            write(iounit, format_9991), path, "Symmetric";
        }
        write(iounit, "(' Matrix types:')");
        if (sord) {
            write(iounit, format_9972);
        } else {
            write(iounit, format_9971);
        }
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9959), 4;
        write(iounit, format_9958), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9957), 7;
        write(iounit, format_9955), 8;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HA")) {
        //
        //        HA: Hermitian,
        //            with Assen Algorithm
        //
        write(iounit, format_9992), path, "Hermitian";
        //
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9972);
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9960), 4;
        write(iounit, format_9959), 5;
        write(iounit, format_9958), 6;
        write(iounit, format_9956), 7;
        write(iounit, format_9957), 8;
        write(iounit, format_9955), 9;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HE")) {
        //
        //        HE: Hermitian indefinite full,
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //
        write(iounit, format_9992), path, "Hermitian";
        //
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9972);
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9960), 4;
        write(iounit, format_9959), 5;
        write(iounit, format_9958), 6;
        write(iounit, format_9956), 7;
        write(iounit, format_9957), 8;
        write(iounit, format_9955), 9;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HR") || Mlsamen(2, p2, "HR")) {
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
        write(iounit, format_9892), path, "Hermitian";
        //
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9972);
        //
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9927), 3;
        write(iounit, format_9928);
        write(iounit, format_9926), 4;
        write(iounit, format_9928);
        write(iounit, format_9960), 5;
        write(iounit, format_9959), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "HP")) {
        //
        //        HP: Hermitian indefinite packed,
        //            with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3, "E")) {
            write(iounit, format_9992), path, "Hermitian";
        } else {
            write(iounit, format_9991), path, "Hermitian";
        }
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9972);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9953), 1;
        write(iounit, format_9961), 2;
        write(iounit, format_9960), 3;
        write(iounit, format_9959), 4;
        write(iounit, format_9958), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9957), 7;
        write(iounit, format_9955), 8;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "TR") || Mlsamen(2, p2, "TP")) {
        //
        //        TR: Triangular full
        //        TP: Triangular packed
        //
        if (Mlsame(c3, "R")) {
            write(iounit, "(/,1x,a3,':  Triangular matrices')"), path;
            subnam = std::string(&path[0]) + std::string("LATRS");
        } else {
            write(iounit, "(/,1x,a3,':  Triangular packed matrices')"), path;
            subnam = std::string(&path[0]) + std::string("LATPS");
        }
        write(iounit, "(' Matrix types for ',a3,' routines:',/,4x,'1. Diagonal',24x,"
                      "'6. Scaled near overflow',/,4x,'2. Random, CNDNUM = 2',14x,"
                      "'7. Identity',/,4x,'3. Random, CNDNUM = sqrt(0.1/EPS)  ',"
                      "'8. Unit triangular, CNDNUM = 2',/,4x,'4. Random, CNDNUM = 0.1/EPS',8x,"
                      "'9. Unit, CNDNUM = sqrt(0.1/EPS)',/,4x,'5. Scaled near underflow',10x,"
                      "'10. Unit, CNDNUM = 0.1/EPS')"),
            path;

        write(iounit, "(' Special types for testing ',a,':',/,3x,"
                      "'11. Matrix elements are O(1), large right hand side',/,3x,"
                      "'12. First diagonal causes overflow,',' offdiagonal column norms < 1',"
                      "/,3x,'13. First diagonal causes overflow,',"
                      "' offdiagonal column norms > 1',/,3x,"
                      "'14. Growth factor underflows, solution does not overflow',/,3x,"
                      "'15. Small diagonal causes gradual overflow',/,3x,"
                      "'16. One zero diagonal element',/,3x,"
                      "'17. Large offdiagonals cause overflow when adding a column',/,3x,"
                      "'18. Unit triangular with large right hand side')"),
            subnam;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9961), 1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, format_9951), subnam, 8;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "TB")) {
        //
        //        TB: Triangular band
        //
        write(iounit, "(/,1x,a3,':  Triangular band matrices')"), path;
        subnam = std::string(&path[0]) + std::string("LATBS");
        write(iounit, "(' Matrix types for ',a3,' routines:',/,4x,'1. Random, CNDNUM = 2',14x,"
                      "'6. Identity',/,4x,'2. Random, CNDNUM = sqrt(0.1/EPS)  ',"
                      "'7. Unit triangular, CNDNUM = 2',/,4x,'3. Random, CNDNUM = 0.1/EPS',8x,"
                      "'8. Unit, CNDNUM = sqrt(0.1/EPS)',/,4x,'4. Scaled near underflow',11x,"
                      "'9. Unit, CNDNUM = 0.1/EPS',/,4x,'5. Scaled near overflow')"),
            path;
        write(iounit, "(' Special types for testing ',a,':',/,3x,"
                      "'10. Matrix elements are O(1), large right hand side',/,3x,"
                      "'11. First diagonal causes overflow,',' offdiagonal column norms < 1',"
                      "/,3x,'12. First diagonal causes overflow,',"
                      "' offdiagonal column norms > 1',/,3x,"
                      "'13. Growth factor underflows, solution does not overflow',/,3x,"
                      "'14. Small diagonal causes gradual overflow',/,3x,"
                      "'15. One zero diagonal element',/,3x,"
                      "'16. Large offdiagonals cause overflow when adding a column',/,3x,"
                      "'17. Unit triangular with large right hand side')"),
            subnam;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9960), 1;
        write(iounit, format_9959), 2;
        write(iounit, format_9958), 3;
        write(iounit, format_9957), 4;
        write(iounit, format_9956), 5;
        write(iounit, format_9955), 6;
        write(iounit, format_9951), subnam, 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "QR")) {
        //
        //        QR decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "QR";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( R - Q'' * A ) / ( M * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( R - Q'' * A ) / ( M * norm(A) * EPS )       [RFPG]')"), 8;
        write(iounit, format_9946), 2;
        write(iounit, format_9944), 3, "M";
        write(iounit, format_9943), 4, "M";
        write(iounit, format_9942), 5, "M";
        write(iounit, format_9941), 6, "M";
        write(iounit, format_9960), 7;
        write(iounit, "(3x,i2,': diagonal is not non-negative')"), 9;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "LQ")) {
        //
        //        LQ decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "LQ";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( L - A * Q'' ) / ( N * norm(A) * EPS )')"), 1;
        write(iounit, format_9945), 2;
        write(iounit, format_9944), 3, "N";
        write(iounit, format_9943), 4, "N";
        write(iounit, format_9942), 5, "N";
        write(iounit, format_9941), 6, "N";
        write(iounit, format_9960), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "QL")) {
        //
        //        QL decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "QL";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( L - Q'' * A ) / ( M * norm(A) * EPS )')"), 1;
        write(iounit, format_9946), 2;
        write(iounit, format_9944), 3, "M";
        write(iounit, format_9943), 4, "M";
        write(iounit, format_9942), 5, "M";
        write(iounit, format_9941), 6, "M";
        write(iounit, format_9960), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "RQ")) {
        //
        //        RQ decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "RQ";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( R - A * Q'' ) / ( N * norm(A) * EPS )')"), 1;
        write(iounit, format_9945), 2;
        write(iounit, format_9944), 3, "N";
        write(iounit, format_9943), 4, "N";
        write(iounit, format_9942), 5, "N";
        write(iounit, format_9941), 6, "N";
        write(iounit, format_9960), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "QP")) {
        //
        //        QR decomposition with column pivoting
        //
        write(iounit, "(/,1x,a3,':  QR factorization with column pivoting')"), path;
        write(iounit, "(' Matrix types (2-6 have condition 1/EPS):',/,4x,'1. Zero matrix',21x,"
                      "'4. First n/2 columns fixed',/,4x,'2. One small eigenvalue',12x,"
                      "'5. Last n/2 columns fixed',/,4x,'3. Geometric distribution',10x,"
                      "'6. Every second column fixed')");
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9940), 1;
        write(iounit, "(3x,i2,': norm( A*P - Q*R )     / ( M * norm(A) * EPS )')"), 2;
        write(iounit, format_9938), 3;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "TZ")) {
        //
        //        TZ:  Trapezoidal
        //
        write(iounit, "(/,1x,a3,':  RQ factorization of trapezoidal matrix')"), path;
        write(iounit, "(' Matrix types (2-3 have condition 1/EPS):',/,4x,'1. Zero matrix',/,"
                      "4x,'2. One small eigenvalue',/,4x,'3. Geometric distribution')");
        write(iounit, "(' Test ratios (1-3: ',a1,'TZRZF):')"), c1;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9940), 1;
        write(iounit, "(3x,i2,': norm( A - R*Q )       / ( M * norm(A) * EPS )')"), 2;
        write(iounit, format_9938), 3;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "LS")) {
        //
        //        LS:  Least Squares driver routines for
        //             LS, LSD, LSS, LSX and LSY.
        //
        write(iounit, "(/,1x,a3,':  Least squares driver routines')"), path;
        write(iounit, "(' Matrix types (1-3: full rank, 4-6: rank deficient):',/,4x,"
                      "'1 and 4. Normal scaling',/,4x,'2 and 5. Scaled near overflow',/,4x,"
                      "'3 and 6. Scaled near underflow')");
        write(iounit, "(' Test ratios:',/,'    (1-2: ',a1,'GELS, 3-6: ',a1,'GELSY, 7-10: ',a1,"
                      "'GELSS, 11-14: ',a1,'GELSD, 15-16: ',a1,'GETSLS)')"),
            c1, c1, c1, c1;
        write(iounit, format_9935), 1;
        write(iounit, "(3x,i2,': norm( (A*X-B)'' *A ) / ',"
                      "'( max(M,N,NRHS) * norm(A) * norm(B) * EPS )',/,7x,"
                      "'if TRANS=''N'' and M.GE.N or TRANS=''T'' and M.LT.N, ','otherwise',/,"
                      "7x,'check if X is in the row space of A or A'' ',"
                      "'(overdetermined case)')"),
            2;
        write(iounit, "(3x,i2,': norm(svd(A)-svd(R)) / ','( min(M,N) * norm(svd(R)) * EPS )')"), 3;
        write(iounit, format_9935), 4;
        write(iounit, "(3x,i2,': norm( (A*X-B)'' *A ) / ',"
                      "'( max(M,N,NRHS) * norm(A) * norm(B) * EPS )')"),
            5;
        write(iounit, "(3x,i2,': Check if X is in the row space of A or A''')"), 6;
        write(iounit, "(3x,' 7-10: same as 3-6',3x,' 11-14: same as 3-6')");
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "LU")) {
        //
        //        LU factorization variants
        //
        write(iounit, "(/,1x,a3,':  LU factorization variants')"), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9979);
        write(iounit, "(' Test ratio:')");
        write(iounit, format_9962), 1;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "CH")) {
        //
        //        Cholesky factorization variants
        //
        write(iounit, "(/,1x,a3,':  Cholesky factorization variants')"), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, "(4x,'1. Diagonal',24x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                      "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = 0.1/EPS',/,3x,"
                      "'*3. First row and column zero',7x,'8. Scaled near underflow',/,3x,"
                      "'*4. Last row and column zero',8x,'9. Scaled near overflow',/,3x,"
                      "'*5. Middle row and column zero',/,3x,"
                      "'(* - tests error exits, no test ratios are computed)')");
        write(iounit, "(' Test ratio:')");
        write(iounit, format_9954), 1;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2, "QS")) {
        //
        //        QR factorization variants
        //
        write(iounit, "(/,1x,a3,':  QR factorization variants')"), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        //
    } else if (Mlsamen(2, p2, "QT")) {
        //
        //        QRT (general matrices)
        //
        write(iounit, "(/,1x,a3,':  QRT factorization for general matrices')"), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( R - Q''*A ) / ( M * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( I - Q''*Q ) / ( M * EPS )')"), 2;
        write(iounit, "(3x,i2,': norm( Q*C - Q*C ) / ( M * norm(C) * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( Q''*C - Q''*C ) / ( M * norm(C) * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( C*Q - C*Q ) / ( M * norm(C) * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( M * norm(C) * EPS )')"), 6;
        //
    } else if (Mlsamen(2, p2, "QX")) {
        //
        //        QRT (triangular-pentagonal)
        //
        write(iounit, "(/,1x,a3,':  QRT factorization for ','triangular-pentagonal matrices')"), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( R - Q''*A ) / ( (M+N) * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( I - Q''*Q ) / ( (M+N) * EPS )')"), 2;
        write(iounit, "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')"), 6;
        //
    } else if (Mlsamen(2, p2, "TQ")) {
        //
        //        QRT (triangular-pentagonal)
        //
        write(iounit, "(/,1x,a3,':  LQT factorization for general matrices')"), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( L - A*Q'' ) / ( (M+N) * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( I - Q*Q'' ) / ( (M+N) * EPS )')"), 2;
        write(iounit, "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')"), 6;
        //
    } else if (Mlsamen(2, p2, "XQ")) {
        //
        //        QRT (triangular-pentagonal)
        //
        write(iounit, "(/,1x,a3,':  LQT factorization for ','triangular-pentagonal matrices')"), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( L - A*Q'' ) / ( (M+N) * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( I - Q*Q'' ) / ( (M+N) * EPS )')"), 2;
        write(iounit, "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')"), 6;
        //
    } else if (Mlsamen(2, p2, "TS")) {
        //
        //        TS:  QR routines for tall-skinny and short-wide matrices
        //
        write(iounit, "(/,1x,a3,':  TS factorization for ',"
                      "'tall-skinny or short-wide matrices')"),
            path;
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( R - Q''*A ) / ( (M+N) * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( I - Q''*Q ) / ( (M+N) * EPS )')"), 2;
        write(iounit, "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')"), 6;
        //
    } else if (Mlsamen(2, p2, "HH")) {
        //
        //        HH:  Householder reconstruction for tall-skinny matrices
        //
        write(iounit, "(/,1x,a3,':  Householder recostruction from TSQR',"
                      "' factorization output ',/,' for tall-skinny matrices.')"),
            path;
        write(iounit, "(' Test ratios:')");
        write(iounit, "(3x,i2,': norm( R - Q''*A ) / ( M * norm(A) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( I - Q''*Q ) / ( M * EPS )')"), 2;
        write(iounit, "(3x,i2,': norm( Q*C - Q*C ) / ( M * norm(C) * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( Q''*C - Q''*C ) / ( M * norm(C) * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( C*Q - C*Q ) / ( M * norm(C) * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( M * norm(C) * EPS )')"), 6;
        //
    } else {
        //
        //        Print error message if no header is available.
        //
        write(iounit, "(/,1x,a3,':  No header available')"), path;
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
    //     CH matrix types
    //
    //     PS matrix types
    //
    //     PB matrix types
    //
    //     SSY, SSR, SSP, CHE, CHR, CHP matrix types
    //
    //     CSY, CSR, CSP matrix types
    //
    //     QR matrix types
    //
    //     QP matrix types
    //
    //     TZ matrix types
    //
    //     LS matrix types
    //
    //     TR, TP matrix types
    //
    //     TB matrix types
    //
    //     Test ratios
    //
    //     End of Alahd
    //
}
