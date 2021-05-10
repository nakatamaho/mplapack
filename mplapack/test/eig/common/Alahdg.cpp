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

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Alahdg(INTEGER const iounit, const char *path) {
    common cmn;
    common_write write(cmn);
    static const char *format_9932 = "(3x,i2,': norm( I - Q''*Q )   / ( N * EPS )')";
    static const char *format_9933 = "(3x,i2,': norm( I - Z''*Z )   / ( P * EPS )')";
    static const char *format_9950 = "(3x,i2,': A-diagonal matrix  B-upper triangular')";
    static const char *format_9951 = "(3x,i2,': A-diagonal matrix  B-lower triangular')";
    static const char *format_9952 = "(3x,i2,': A-upper triangular B-upper triangular')";
    static const char *format_9953 = "(3x,i2,': A-lower triangular B-diagonal triangular')";
    static const char *format_9954 = "(3x,i2,': A-lower triangular B-upper triangular')";
    static const char *format_9955 = "(3x,i2,': Random matrices cond(A)=100, cond(B)=10,')";
    static const char *format_9956 = "(3x,i2,': Random matrices cond(A)= sqrt( 0.1/EPS ) ',"
                                     "'cond(B)= sqrt( 0.1/EPS )')";
    static const char *format_9957 = "(3x,i2,': Random matrices cond(A)= 0.1/EPS ','cond(B)= 0.1/EPS')";
    static const char *format_9961 = "(3x,i2,': Matrix scaled near underflow limit')";
    static const char *format_9962 = "(3x,i2,': Matrix scaled near overflow limit')";
    static const char *format_9999 = "(1x,a)";
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
    char c2[3];
    c2[0] = path[0];
    c2[1] = path[1];
    c2[2] = path[2];
    //
    //     First line describing matrices in this path
    //
    INTEGER itype = 0;
    if (Mlsamen(3, c2, "GQR")) {
        itype = 1;
        write(iounit, "(/,1x,a3,': GQR factorization of general matrices')"), path;
    } else if (Mlsamen(3, c2, "GRQ")) {
        itype = 2;
        write(iounit, "(/,1x,a3,': GRQ factorization of general matrices')"), path;
    } else if (Mlsamen(3, c2, "LSE")) {
        itype = 3;
        write(iounit, "(/,1x,a3,': LSE Problem')"), path;
    } else if (Mlsamen(3, c2, "GLM")) {
        itype = 4;
        write(iounit, "(/,1x,a3,': GLM Problem')"), path;
    } else if (Mlsamen(3, c2, "GSV")) {
        itype = 5;
        write(iounit, "(/,1x,a3,': Generalized Singular Value Decomposition')"), path;
    } else if (Mlsamen(3, c2, "CSD")) {
        itype = 6;
        write(iounit, "(/,1x,a3,': CS Decomposition')"), path;
    }
    //
    //     Matrix types
    //
    write(iounit, format_9999), "Matrix types: ";
    //
    if (itype == 1) {
        write(iounit, format_9950), 1;
        write(iounit, format_9952), 2;
        write(iounit, format_9954), 3;
        write(iounit, format_9955), 4;
        write(iounit, format_9956), 5;
        write(iounit, format_9957), 6;
        write(iounit, format_9961), 7;
        write(iounit, format_9962), 8;
    } else if (itype == 2) {
        write(iounit, format_9951), 1;
        write(iounit, format_9953), 2;
        write(iounit, format_9954), 3;
        write(iounit, format_9955), 4;
        write(iounit, format_9956), 5;
        write(iounit, format_9957), 6;
        write(iounit, format_9961), 7;
        write(iounit, format_9962), 8;
    } else if (itype == 3) {
        write(iounit, format_9950), 1;
        write(iounit, format_9952), 2;
        write(iounit, format_9954), 3;
        write(iounit, format_9955), 4;
        write(iounit, format_9955), 5;
        write(iounit, format_9955), 6;
        write(iounit, format_9955), 7;
        write(iounit, format_9955), 8;
    } else if (itype == 4) {
        write(iounit, format_9951), 1;
        write(iounit, format_9953), 2;
        write(iounit, format_9954), 3;
        write(iounit, format_9955), 4;
        write(iounit, format_9955), 5;
        write(iounit, format_9955), 6;
        write(iounit, format_9955), 7;
        write(iounit, format_9955), 8;
    } else if (itype == 5) {
        write(iounit, format_9950), 1;
        write(iounit, format_9952), 2;
        write(iounit, format_9954), 3;
        write(iounit, format_9955), 4;
        write(iounit, format_9956), 5;
        write(iounit, format_9957), 6;
        write(iounit, "(3x,i2,': Random matrices cond(A)= sqrt( 0.1/EPS ) ',"
                      "'cond(B)=  0.1/EPS ')"),
            7;
        write(iounit, "(3x,i2,': Random matrices cond(A)= 0.1/EPS ',"
                      "'cond(B)=  sqrt( 0.1/EPS )')"),
            8;
    } else if (itype == 6) {
        write(iounit, "(3x,i2,': Random orthogonal matrix (Haar measure)')"), 1;
        write(iounit, "(3x,i2,': Nearly orthogonal matrix with uniformly ',"
                      "'distributed angles atan2( S, C ) in CS decomposition')"),
            2;
        write(iounit, "(3x,i2,': Random orthogonal matrix with clustered ',"
                      "'angles atan2( S, C ) in CS decomposition')"),
            3;
    }
    //
    //     Tests performed
    //
    write(iounit, format_9999), "Test ratios: ";
    //
    if (itype == 1) {
        //
        //        GQR decomposition of rectangular matrices
        //
        write(iounit, "(3x,i2,': norm( R - Q'' * A ) / ( min( N, M )*norm( A )','* EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( T * Z - Q'' * B )  / ( min(P,N)*norm(B)','* EPS )')"), 2;
        write(iounit, format_9932), 3;
        write(iounit, format_9933), 4;
    } else if (itype == 2) {
        //
        //        GRQ decomposition of rectangular matrices
        //
        write(iounit, "(3x,i2,': norm( R - A * Q'' ) / ( min( N,M )*norm(A) * ','EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( T * Q - Z'' * B )  / ( min( P,N ) * nor','m(B)*EPS )')"), 2;
        write(iounit, format_9932), 3;
        write(iounit, format_9933), 4;
    } else if (itype == 3) {
        //
        //        LSE Problem
        //
        write(iounit, "(3x,i2,': norm( A*x - c )  / ( norm(A)*norm(x) * EPS )')"), 1;
        write(iounit, "(3x,i2,': norm( B*x - d )  / ( norm(B)*norm(x) * EPS )')"), 2;
    } else if (itype == 4) {
        //
        //        GLM Problem
        //
        write(iounit, "(3x,i2,': norm( d - A*x - B*y ) / ( (norm(A)+norm(B) )*',"
                      "'(norm(x)+norm(y))*EPS )')"),
            1;
    } else if (itype == 5) {
        //
        //        GSVD
        //
        write(iounit, "(3x,i2,': norm( U'' * A * Q - D1 * R ) / ( min( M, N )*',"
                      "'norm( A ) * EPS )')"),
            1;
        write(iounit, "(3x,i2,': norm( V'' * B * Q - D2 * R ) / ( min( P, N )*',"
                      "'norm( B ) * EPS )')"),
            2;
        write(iounit, "(3x,i2,': norm( I - U''*U )   / ( M * EPS )')"), 3;
        write(iounit, "(3x,i2,': norm( I - V''*V )   / ( P * EPS )')"), 4;
        write(iounit, "(3x,i2,': norm( I - Q''*Q )   / ( N * EPS )')"), 5;
    } else if (itype == 6) {
        //
        //        CSD
        //
        write(iounit, "(3x,'2-by-2 CSD')");
        write(iounit, "(3x,i2,': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)',"
                      "' * max(norm(I-X''*X),EPS) )')"),
            1;
        write(iounit, "(3x,i2,': norm( U1'' * X12 * V2-(-S)) / ( max(  P,',"
                      "'M-Q) * max(norm(I-X''*X),EPS) )')"),
            2;
        write(iounit, "(3x,i2,': norm( U2'' * X21 * V1 - S ) / ( max(M-P,',"
                      "'  Q) * max(norm(I-X''*X),EPS) )')"),
            3;
        write(iounit, "(3x,i2,': norm( U2'' * X22 * V2 - C ) / ( max(M-P,',"
                      "'M-Q) * max(norm(I-X''*X),EPS) )')"),
            4;
        write(iounit, "(3x,i2,': norm( I - U1''*U1 ) / (   P   * EPS )')"), 5;
        write(iounit, "(3x,i2,': norm( I - U2''*U2 ) / ( (M-P) * EPS )')"), 6;
        write(iounit, "(3x,i2,': norm( I - V1''*V1 ) / (   Q   * EPS )')"), 7;
        write(iounit, "(3x,i2,': norm( I - V2''*V2 ) / ( (M-Q) * EPS )')"), 8;
        write(iounit, "(3x,i2,': principal angle ordering ( 0 or ULP )')"), 9;
        write(iounit, "(3x,'2-by-1 CSD')");
        write(iounit, "(3x,i2,': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)',"
                      "' * max(norm(I-X''*X),EPS) )')"),
            10;
        write(iounit, "(3x,i2,': norm( U2'' * X21 * V1 - S ) / ( max(  M-P,',"
                      "'Q) * max(norm(I-X''*X),EPS) )')"),
            11;
        write(iounit, "(3x,i2,': norm( I - U1''*U1 ) / (   P   * EPS )')"), 12;
        write(iounit, "(3x,i2,': norm( I - U2''*U2 ) / ( (M-P) * EPS )')"), 13;
        write(iounit, "(3x,i2,': norm( I - V1''*V1 ) / (   Q   * EPS )')"), 14;
        write(iounit, "(3x,i2,': principal angle ordering ( 0 or ULP )')"), 15;
    }
    //
    //     GQR test ratio
    //
    //     GRQ test ratio
    //
    //     LSE test ratio
    //
    //     GLM test ratio
    //
    //     GSVD test ratio
    //
    //     CSD test ratio
    //
    //     End of Alahdg
    //
}
