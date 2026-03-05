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

// Derived from LAPACK routine ALAHDG.
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
#include <mplapack_eig.h>

void Alahdg(INTEGER const iounit, fem::str_cref path) {
    common cmn;
    common_write write(cmn);
    //
    static const char *format_9999 = "(1x,a)";
    static const char *format_9991 = "(/,1x,a3,': GQR factorization of general matrices')";
    static const char *format_9992 = "(/,1x,a3,': GRQ factorization of general matrices')";
    static const char *format_9993 = "(/,1x,a3,': LSE Problem')";
    static const char *format_9994 = "(/,1x,a3,': GLM Problem')";
    static const char *format_9995 = "(/,1x,a3,': Generalized Singular Value Decomposition')";
    static const char *format_9996 = "(/,1x,a3,': CS Decomposition')";
    //
    static const char *format_9950 = "(3x,i2,': A-diagonal matrix  B-upper triangular')";
    static const char *format_9951 = "(3x,i2,': A-diagonal matrix  B-lower triangular')";
    static const char *format_9952 = "(3x,i2,': A-upper triangular B-upper triangular')";
    static const char *format_9953 = "(3x,i2,': A-lower triangular B-diagonal triangular')";
    static const char *format_9954 = "(3x,i2,': A-lower triangular B-upper triangular')";
    //
    static const char *format_9955 = "(3x,i2,': Random matrices cond(A)=100, cond(B)=10,')";
    //
    static const char *format_9956 = "(3x,i2,': Random matrices cond(A)= sqrt( 0.1/EPS ) ',"
                                     "'cond(B)= sqrt( 0.1/EPS )')";
    static const char *format_9957 = "(3x,i2,': Random matrices cond(A)= 0.1/EPS ','cond(B)= 0.1/EPS')";
    static const char *format_9959 = "(3x,i2,': Random matrices cond(A)= sqrt( 0.1/EPS ) ',"
                                     "'cond(B)=  0.1/EPS ')";
    static const char *format_9960 = "(3x,i2,': Random matrices cond(A)= 0.1/EPS ','cond(B)=  sqrt( 0.1/EPS )')";
    //
    static const char *format_9961 = "(3x,i2,': Matrix scaled near underflow limit')";
    static const char *format_9962 = "(3x,i2,': Matrix scaled near overflow limit')";
    static const char *format_9963 = "(3x,i2,': Random orthogonal matrix (Haar measure)')";
    static const char *format_9964 = "(3x,i2,': Nearly orthogonal matrix with uniformly ',"
                                     "'distributed angles atan2( S, C ) in CS decomposition')";
    static const char *format_9965 = "(3x,i2,': Random orthogonal matrix with clustered ',"
                                     "'angles atan2( S, C ) in CS decomposition')";
    //
    //
    // GQR test ratio
    //
    static const char *format_9930 = "(3x,i2,': norm( R - Q'' * A ) / ( min( N, M )*norm( A )','* EPS )')";
    static const char *format_9931 = "(3x,i2,': norm( T * Z - Q'' * B )  / ( min(P,N)*norm(B)','* EPS )')";
    static const char *format_9932 = "(3x,i2,': norm( I - Q''*Q )   / ( N * EPS )')";
    static const char *format_9933 = "(3x,i2,': norm( I - Z''*Z )   / ( P * EPS )')";
    //
    // GRQ test ratio
    //
    static const char *format_9934 = "(3x,i2,': norm( R - A * Q'' ) / ( min( N,M )*norm(A) * ','EPS )')";
    static const char *format_9935 = "(3x,i2,': norm( T * Q - Z'' * B )  / ( min( P,N ) * nor','m(B)*EPS )')";
    //
    // LSE test ratio
    //
    static const char *format_9937 = "(3x,i2,': norm( A*x - c )  / ( norm(A)*norm(x) * EPS )')";
    static const char *format_9938 = "(3x,i2,': norm( B*x - d )  / ( norm(B)*norm(x) * EPS )')";
    //
    // GLM test ratio
    //
    static const char *format_9939 = "(3x,i2,': norm( d - A*x - B*y ) / ( (norm(A)+norm(B) )*',"
                                     "'(norm(x)+norm(y))*EPS )')";
    //
    // GSVD test ratio
    //
    static const char *format_9940 = "(3x,i2,': norm( U'' * A * Q - D1 * R ) / ( min( M, N )*',"
                                     "'norm( A ) * EPS )')";
    static const char *format_9941 = "(3x,i2,': norm( V'' * B * Q - D2 * R ) / ( min( P, N )*',"
                                     "'norm( B ) * EPS )')";
    static const char *format_9942 = "(3x,i2,': norm( I - U''*U )   / ( M * EPS )')";
    static const char *format_9943 = "(3x,i2,': norm( I - V''*V )   / ( P * EPS )')";
    static const char *format_9944 = "(3x,i2,': norm( I - Q''*Q )   / ( N * EPS )')";
    //
    // CSD test ratio
    //
    static const char *format_9910 = "(3x,'2-by-2 CSD')";
    static const char *format_9911 = "(3x,i2,': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)',"
                                     "' * max(norm(I-X''*X),EPS) )')";
    static const char *format_9912 = "(3x,i2,': norm( U1'' * X12 * V2-(-S)) / ( max(  P,',"
                                     "'M-Q) * max(norm(I-X''*X),EPS) )')";
    static const char *format_9913 = "(3x,i2,': norm( U2'' * X21 * V1 - S ) / ( max(M-P,',"
                                     "'  Q) * max(norm(I-X''*X),EPS) )')";
    static const char *format_9914 = "(3x,i2,': norm( U2'' * X22 * V2 - C ) / ( max(M-P,',"
                                     "'M-Q) * max(norm(I-X''*X),EPS) )')";
    static const char *format_9915 = "(3x,i2,': norm( I - U1''*U1 ) / (   P   * EPS )')";
    static const char *format_9916 = "(3x,i2,': norm( I - U2''*U2 ) / ( (M-P) * EPS )')";
    static const char *format_9917 = "(3x,i2,': norm( I - V1''*V1 ) / (   Q   * EPS )')";
    static const char *format_9918 = "(3x,i2,': norm( I - V2''*V2 ) / ( (M-Q) * EPS )')";
    static const char *format_9919 = "(3x,i2,': principal angle ordering ( 0 or ULP )')";
    static const char *format_9920 = "(3x,'2-by-1 CSD')";
    static const char *format_9921 = "(3x,i2,': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)',"
                                     "' * max(norm(I-X''*X),EPS) )')";
    static const char *format_9922 = "(3x,i2,': norm( U2'' * X21 * V1 - S ) / ( max(  M-P,',"
                                     "'Q) * max(norm(I-X''*X),EPS) )')";
    static const char *format_9923 = "(3x,i2,': norm( I - U1''*U1 ) / (   P   * EPS )')";
    static const char *format_9924 = "(3x,i2,': norm( I - U2''*U2 ) / ( (M-P) * EPS )')";
    static const char *format_9925 = "(3x,i2,': norm( I - V1''*V1 ) / (   Q   * EPS )')";
    static const char *format_9926 = "(3x,i2,': principal angle ordering ( 0 or ULP )')";
    //
    if (iounit <= 0) {
        return;
    }
    fem::str<3> c2 = path(1, 3);
    //
    // First line describing matrices in this path
    //
    INTEGER itype = 0;
    if (Mlsamen(3, c2.elems, "GQR")) {
        itype = 1;
        write(iounit, format_9991), path;
    } else if (Mlsamen(3, c2.elems, "GRQ")) {
        itype = 2;
        write(iounit, format_9992), path;
    } else if (Mlsamen(3, c2.elems, "LSE")) {
        itype = 3;
        write(iounit, format_9993), path;
    } else if (Mlsamen(3, c2.elems, "GLM")) {
        itype = 4;
        write(iounit, format_9994), path;
    } else if (Mlsamen(3, c2.elems, "GSV")) {
        itype = 5;
        write(iounit, format_9995), path;
    } else if (Mlsamen(3, c2.elems, "CSD")) {
        itype = 6;
        write(iounit, format_9996), path;
    }
    //
    // Matrix types
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
        write(iounit, format_9959), 7;
        write(iounit, format_9960), 8;
    } else if (itype == 6) {
        write(iounit, format_9963), 1;
        write(iounit, format_9964), 2;
        write(iounit, format_9965), 3;
    }
    //
    // Tests performed
    //
    write(iounit, format_9999), "Test ratios: ";
    //
    if (itype == 1) {
        //
        // GQR decomposition of rectangular matrices
        //
        write(iounit, format_9930), 1;
        write(iounit, format_9931), 2;
        write(iounit, format_9932), 3;
        write(iounit, format_9933), 4;
    } else if (itype == 2) {
        //
        // GRQ decomposition of rectangular matrices
        //
        write(iounit, format_9934), 1;
        write(iounit, format_9935), 2;
        write(iounit, format_9932), 3;
        write(iounit, format_9933), 4;
    } else if (itype == 3) {
        //
        // LSE Problem
        //
        write(iounit, format_9937), 1;
        write(iounit, format_9938), 2;
    } else if (itype == 4) {
        //
        // GLM Problem
        //
        write(iounit, format_9939), 1;
    } else if (itype == 5) {
        //
        // GSVD
        //
        write(iounit, format_9940), 1;
        write(iounit, format_9941), 2;
        write(iounit, format_9942), 3;
        write(iounit, format_9943), 4;
        write(iounit, format_9944), 5;
    } else if (itype == 6) {
        //
        // CSD
        //
        write(iounit, format_9910);
        write(iounit, format_9911), 1;
        write(iounit, format_9912), 2;
        write(iounit, format_9913), 3;
        write(iounit, format_9914), 4;
        write(iounit, format_9915), 5;
        write(iounit, format_9916), 6;
        write(iounit, format_9917), 7;
        write(iounit, format_9918), 8;
        write(iounit, format_9919), 9;
        write(iounit, format_9920);
        write(iounit, format_9921), 10;
        write(iounit, format_9922), 11;
        write(iounit, format_9923), 12;
        write(iounit, format_9924), 13;
        write(iounit, format_9925), 14;
        write(iounit, format_9926), 15;
    }
    //
    // End of Alahdg
    //
}
