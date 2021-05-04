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

void Rlahd2(INTEGER const iounit, const char *path) {
    common_write write(cmn);
    static const char *format_9969 = "(/,' Test ratios:  ','(B: upper bidiagonal, Q and P: ',a10,/,16x,"
                                     "'C: m x nrhs, PT = P'', Y = Q'' C)',/,"
                                     "' 1: norm( A - Q B PT ) / ( norm(A) max(m,n) ulp )',/,"
                                     "' 2: norm( I - Q'' Q )   / ( m ulp )',/,"
                                     "' 3: norm( I - PT PT'' )   / ( n ulp )',/,"
                                     "' 4: norm( Y - Q'' C )   / ( norm(Y) max(m,nrhs) ulp )')";
    static const char *format_9970 = "(' Matrix types (see xCHKBB for details):',/,' Diagonal matrices:',/,"
                                     "'   1: Zero',28x,' 5: Clustered entries',/,'   2: Identity',24x,"
                                     "' 6: Large, evenly spaced entries',/,'   3: Evenly spaced entries',11x,"
                                     "' 7: Small, evenly spaced entries',/,"
                                     "'   4: Geometrically spaced entries',/,' General matrices:',/,"
                                     "'   8: Evenly spaced sing. vals.',7x,"
                                     "'12: Small, evenly spaced sing vals',/,"
                                     "'   9: Geometrically spaced sing vals  ','13: Random, O(1) entries',/,"
                                     "'  10: Clustered sing. vals.',11x,'14: Random, scaled near overflow',/,"
                                     "'  11: Large, evenly spaced sing vals  ',"
                                     "'15: Random, scaled near underflow')";
    static const char *format_9971 = "('   1: norm( A - Q B P'' ) / ( norm(A) max(m,n) ulp )',/,"
                                     "'   2: norm( I - Q'' Q )   / ( m ulp )',/,"
                                     "'   3: norm( I - P'' P )   / ( n ulp )',/,"
                                     "'   4: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )',/,"
                                     "'   5: norm( Y - U Z )    / ','( norm(Z) max(min(m,n),k) ulp )',/,"
                                     "'   6: norm( I - U'' U )   / ( min(m,n) ulp )',/,"
                                     "'   7: norm( I - V'' V )   / ( min(m,n) ulp )',/,"
                                     "'   8: Test ordering of S  (0 if nondecreasing, 1/ulp ',' otherwise)',/,"
                                     "'   9: norm( S - S1 )     / ( norm(S) ulp ),',' where S1 is computed',/,"
                                     "43x,' without computing U and V''',/,'  10: Sturm sequence test ',"
                                     "'(0 if sing. vals of B within THRESH of S)',/,"
                                     "'  11: norm( A - (QU) S (V'' P'') ) / ','( norm(A) max(m,n) ulp )',/,"
                                     "'  12: norm( X - (QU) Z )         / ( |X| max(M,k) ulp )',/,"
                                     "'  13: norm( I - (QU)''(QU) )      / ( M ulp )',/,"
                                     "'  14: norm( I - (V'' P'') (P V) )  / ( N ulp )',/,"
                                     "'  15: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )',/,"
                                     "'  16: norm( I - U'' U )   / ( min(m,n) ulp )',/,"
                                     "'  17: norm( I - V'' V )   / ( min(m,n) ulp )',/,"
                                     "'  18: Test ordering of S  (0 if nondecreasing, 1/ulp ',' otherwise)',/,"
                                     "'  19: norm( S - S1 )     / ( norm(S) ulp ),',' where S1 is computed',/,"
                                     "43x,' without computing U and V''',/,"
                                     "'  20: norm( B - U S V'' )  / ( norm(B) min(m,n) ulp )','  DBDSVX(V,A)',"
                                     "/,'  21: norm( I - U'' U )    / ( min(m,n) ulp )',/,"
                                     "'  22: norm( I - V'' V )    / ( min(m,n) ulp )',/,"
                                     "'  23: Test ordering of S  (0 if nondecreasing, 1/ulp ',' otherwise)',/,"
                                     "'  24: norm( S - S1 )      / ( norm(S) ulp ),',' where S1 is computed',/,"
                                     "44x,' without computing U and V''',/,"
                                     "'  25: norm( S - U'' B V ) / ( norm(B) n ulp )','  DBDSVX(V,I)',/,"
                                     "'  26: norm( I - U'' U )    / ( min(m,n) ulp )',/,"
                                     "'  27: norm( I - V'' V )    / ( min(m,n) ulp )',/,"
                                     "'  28: Test ordering of S  (0 if nondecreasing, 1/ulp ',' otherwise)',/,"
                                     "'  29: norm( S - S1 )      / ( norm(S) ulp ),',' where S1 is computed',/,"
                                     "44x,' without computing U and V''',/,"
                                     "'  30: norm( S - U'' B V ) / ( norm(B) n ulp )','  DBDSVX(V,V)',/,"
                                     "'  31: norm( I - U'' U )    / ( min(m,n) ulp )',/,"
                                     "'  32: norm( I - V'' V )    / ( min(m,n) ulp )',/,"
                                     "'  33: Test ordering of S  (0 if nondecreasing, 1/ulp ',' otherwise)',/,"
                                     "'  34: norm( S - S1 )      / ( norm(S) ulp ),',' where S1 is computed',/,"
                                     "44x,' without computing U and V''')";
    static const char *format_9972 = "(/,' Test ratios:  ','(B: bidiagonal, S: diagonal, Q, P, U, and V: ',a10,"
                                     "/,16x,'X: m x nrhs, Y = Q'' X, and Z = U'' Y)')";
    static const char *format_9973 = "(' Matrix types (see xCHKBD for details):',/,' Diagonal matrices:',/,"
                                     "'   1: Zero',28x,' 5: Clustered entries',/,'   2: Identity',24x,"
                                     "' 6: Large, evenly spaced entries',/,'   3: Evenly spaced entries',11x,"
                                     "' 7: Small, evenly spaced entries',/,"
                                     "'   4: Geometrically spaced entries',/,' General matrices:',/,"
                                     "'   8: Evenly spaced sing. vals.',7x,"
                                     "'12: Small, evenly spaced sing vals',/,"
                                     "'   9: Geometrically spaced sing vals  ','13: Random, O(1) entries',/,"
                                     "'  10: Clustered sing. vals.',11x,'14: Random, scaled near overflow',/,"
                                     "'  11: Large, evenly spaced sing vals  ',"
                                     "'15: Random, scaled near underflow')";
    static const char *format_9978 = "(' Dense or Banded ',a,' Matrices: ',/,"
                                     "'  8=Evenly spaced eigenvals.         ',"
                                     "' 15=Matrix with small random entries.',/,"
                                     "'  9=Geometrically spaced eigenvals.  ',"
                                     "' 16=Evenly spaced eigenvals, KA=1, KB=1.',/,"
                                     "' 10=Clustered eigenvalues.           ',"
                                     "' 17=Evenly spaced eigenvals, KA=2, KB=1.',/,"
                                     "' 11=Large, evenly spaced eigenvals.  ',"
                                     "' 18=Evenly spaced eigenvals, KA=2, KB=2.',/,"
                                     "' 12=Small, evenly spaced eigenvals.  ',"
                                     "' 19=Evenly spaced eigenvals, KA=3, KB=1.',/,"
                                     "' 13=Matrix with random O(1) entries. ',"
                                     "' 20=Evenly spaced eigenvals, KA=3, KB=2.',/,"
                                     "' 14=Matrix with large random entries.',"
                                     "' 21=Evenly spaced eigenvals, KA=3, KB=3.')";
    static const char *format_9979 = "(/,' Special Matrices:',/,'  1=Zero matrix.             ','           ',"
                                     "'  5=Diagonal: clustered entries.',/,'  2=',"
                                     "'Identity matrix.                    ','  6=Diagonal: lar',"
                                     "'ge, evenly spaced.',/,'  3=Diagonal: evenly spaced entri','es.    ',"
                                     "'  7=Diagonal: small, evenly spaced.',/,'  4=D',"
                                     "'iagonal: geometr. spaced entries.')";
    static const char *format_9980 = "(' Matrix types (see xDRVSG for details): ')";
    static const char *format_9981 = "(' Dense ',a,' Matrices:',/,'  8=Evenly spaced eigen',"
                                     "'vals.            ',' 12=Small, evenly spaced eigenvals.',/,"
                                     "'  9=Geometrically spaced eigenvals.     ',' 13=Matrix ',"
                                     "'with random O(1) entries.',/,' 10=Clustered eigenvalues.',"
                                     "'              ',' 14=Matrix with large random entries.',/,"
                                     "' 11=Large, evenly spaced eigenvals.     ',' 15=Matrix ',"
                                     "'with small random entries.')";
    static const char *format_9982 = "(/,' Special Matrices:',/,'  1=Zero matrix.             ','           ',"
                                     "'  5=Diagonal: clustered entries.',/,'  2=',"
                                     "'Identity matrix.                    ','  6=Diagonal: lar',"
                                     "'ge, evenly spaced.',/,'  3=Diagonal: evenly spaced entri','es.    ',"
                                     "'  7=Diagonal: small, evenly spaced.',/,'  4=D',"
                                     "'iagonal: geometr. spaced entries.')";
    static const char *format_9983 = "(' Matrix types (see xDRVST for details): ')";
    static const char *format_9984 = "(/,' Tests performed:   ','(H is Hessenberg, T is Schur,',"
                                     "' U and Z are ',a,',',/,20x,a,', W is a diagonal matr',"
                                     "'ix of eigenvalues,',/,20x,'L and R are the left and rig',"
                                     "'ht eigenvector matrices)',/,'  1 = | A - U H U',a1,' |',"
                                     "' / ( |A| n ulp )         ','  2 = | I - U U',a1,' | / ','( n ulp )',/,"
                                     "'  3 = | H - Z T Z',a1,' | / ( |H| n ulp ',')         ',"
                                     "'  4 = | I - Z Z',a1,' | / ( n ulp )',/,'  5 = | A - UZ T (UZ)',a1,"
                                     "' | / ( |A| n ulp )     ','  6 = | I - UZ (UZ)',a1,' | / ( n ulp )',/,"
                                     "'  7 = | T(','e.vects.) - T(no e.vects.) | / ( |T| ulp )',/,'  8 = | W',"
                                     "'(e.vects.) - W(no e.vects.) | / ( |W| ulp )',/,'  9 = | ',"
                                     "'TR - RW | / ( |T| |R| ulp )     ',' 10 = | LT - WL | / (',"
                                     "' |T| |L| ulp )',/,' 11= |HX - XW| / (|H| |X| ulp)  (inv.','it)',"
                                     "' 12= |YH - WY| / (|H| |Y| ulp)  (inv.it)')";
    static const char *format_9985 = "(' 19=Matrix with random O(1) entries.    ',' 21=Matrix ',"
                                     "'with small random entries.',/,' 20=Matrix with large ran',"
                                     "'dom entries.   ')";
    static const char *format_9986 = "(' Dense, Non-Symmetric Matrices:',/,'  9=Well-cond., ev',"
                                     "'enly spaced eigenvals.',' 14=Ill-cond., geomet. spaced e','igenals.',/,"
                                     "' 10=Well-cond., geom. spaced eigenvals. ',"
                                     "' 15=Ill-conditioned, clustered e.vals.',/,' 11=Well-cond',"
                                     "'itioned, clustered e.vals. ',' 16=Ill-cond., random comp','lex ',a6,/,"
                                     "' 12=Well-cond., random complex ',a6,'   ',"
                                     "' 17=Ill-cond., large rand. complx ',a4,/,' 13=Ill-condi',"
                                     "'tioned, evenly spaced.     ',' 18=Ill-cond., small rand.',' complx ',a4)";
    static const char *format_9987 = "(/,' Special Matrices:',/,'  1=Zero matrix.             ','           ',"
                                     "'  5=Diagonal: geometr. spaced entries.',/,"
                                     "'  2=Identity matrix.                    ','  6=Diagona',"
                                     "'l: clustered entries.',/,'  3=Transposed Jordan block.  ','          ',"
                                     "'  7=Diagonal: large, evenly spaced.',/,'  ',"
                                     "'4=Diagonal: evenly spaced entries.    ','  8=Diagonal: s',"
                                     "'mall, evenly spaced.')";
    static const char *format_9988 = "(' Matrix types (see xCHKHS for details): ')";
    static const char *format_9999 = "(1x,a3,':  no header available')";
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
    bool sord = Mlsame(path, "S") || Mlsame(path, "D");
    bool corz = Mlsame(path, "C") || Mlsame(path, "Z");
    if (!sord && !corz) {
        write(iounit, format_9999), path;
    }
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
    //
    INTEGER j = 0;
    if (Mlsamen(2, c2, "HS")) {
        if (sord) {
            //
            //           Real Non-symmetric Eigenvalue Problem:
            //
            write(iounit, "(/,1x,a3,' -- Real Non-symmetric eigenvalue problem')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9988);
            write(iounit, format_9987);
            write(iounit, format_9986), "pairs ", "pairs ", "prs.", "prs.";
            write(iounit, format_9985);
            //
            //           Tests performed
            //
            {
                write_loop wloop(cmn, iounit, format_9984);
                wloop, "orthogonal", "'=transpose";
                for (j = 1; j <= 6; j = j + 1) {
                    wloop, "'";
                }
            }
            //
        } else {
            //
            //           Complex Non-symmetric Eigenvalue Problem:
            //
            write(iounit, "(/,1x,a3,' -- Complex Non-symmetric eigenvalue problem')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9988);
            write(iounit, format_9987);
            write(iounit, format_9986), "e.vals", "e.vals", "e.vs", "e.vs";
            write(iounit, format_9985);
            //
            //           Tests performed
            //
            {
                write_loop wloop(cmn, iounit, format_9984);
                wloop, "unitary", "*=conj.transp.";
                for (j = 1; j <= 6; j = j + 1) {
                    wloop, "*";
                }
            }
        }
        //
    } else if (Mlsamen(2, c2, "ST")) {
        //
        if (sord) {
            //
            //           Real Symmetric Eigenvalue Problem:
            //
            write(iounit, "(/,1x,a3,' -- Real Symmetric eigenvalue problem')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9983);
            write(iounit, format_9982);
            write(iounit, format_9981), "Symmetric";
            //
            //           Tests performed
            //
            write(iounit, "(/,' Tests performed:  See sdrvst.f')");
            //
        } else {
            //
            //           Complex Hermitian Eigenvalue Problem:
            //
            write(iounit, "(/,1x,a3,' -- Complex Hermitian eigenvalue problem')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9983);
            write(iounit, format_9982);
            write(iounit, format_9981), "Hermitian";
            //
            //           Tests performed
            //
            write(iounit, "(/,' Tests performed:  See cdrvst.f')");
        }
        //
    } else if (Mlsamen(2, c2, "SG")) {
        //
        if (sord) {
            //
            //           Real Symmetric Generalized Eigenvalue Problem:
            //
            write(iounit, "(/,1x,a3,' -- Real Symmetric Generalized eigenvalue ','problem')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9980);
            write(iounit, format_9979);
            write(iounit, format_9978), "Symmetric";
            //
            //           Tests performed
            //
            write(iounit, "(/,' Tests performed:   ',/,"
                          "'( For each pair (A,B), where A is of the given type ',/,"
                          "' and B is a random well-conditioned matrix. D is ',/,"
                          "' diagonal, and Z is orthogonal. )',/,"
                          "' 1 = Rsygv, with ITYPE=1 and UPLO=''U'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 2 = Rspgv, with ITYPE=1 and UPLO=''U'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 3 = Rsbgv, with ITYPE=1 and UPLO=''U'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 4 = Rsygv, with ITYPE=1 and UPLO=''L'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 5 = Rspgv, with ITYPE=1 and UPLO=''L'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 6 = Rsbgv, with ITYPE=1 and UPLO=''L'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ')");
            write(iounit, "(' 7 = Rsygv, with ITYPE=2 and UPLO=''U'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 8 = Rspgv, with ITYPE=2 and UPLO=''U'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 9 = Rspgv, with ITYPE=2 and UPLO=''L'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'10 = Rspgv, with ITYPE=2 and UPLO=''L'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'11 = Rsygv, with ITYPE=3 and UPLO=''U'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'12 = Rspgv, with ITYPE=3 and UPLO=''U'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'13 = Rsygv, with ITYPE=3 and UPLO=''L'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'14 = Rspgv, with ITYPE=3 and UPLO=''L'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ')");
            //
        } else {
            //
            //           Complex Hermitian Generalized Eigenvalue Problem:
            //
            write(iounit, "(/,1x,a3,' -- Complex Hermitian Generalized eigenvalue ','problem')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9980);
            write(iounit, format_9979);
            write(iounit, format_9978), "Hermitian";
            //
            //           Tests performed
            //
            write(iounit, "(/,' Tests performed:   ',/,"
                          "'( For each pair (A,B), where A is of the given type ',/,"
                          "' and B is a random well-conditioned matrix. D is ',/,"
                          "' diagonal, and Z is unitary. )',/,"
                          "' 1 = Chegv, with ITYPE=1 and UPLO=''U'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 2 = Chpgv, with ITYPE=1 and UPLO=''U'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 3 = Chbgv, with ITYPE=1 and UPLO=''U'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 4 = Chegv, with ITYPE=1 and UPLO=''L'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 5 = Chpgv, with ITYPE=1 and UPLO=''L'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 6 = Chbgv, with ITYPE=1 and UPLO=''L'':',"
                          "'  | A Z - B Z D | / ( |A| |Z| n ulp )     ')");
            write(iounit, "(' 7 = Chegv, with ITYPE=2 and UPLO=''U'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 8 = Chpgv, with ITYPE=2 and UPLO=''U'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "' 9 = Chpgv, with ITYPE=2 and UPLO=''L'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'10 = Chpgv, with ITYPE=2 and UPLO=''L'':',"
                          "'  | A B Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'11 = Chegv, with ITYPE=3 and UPLO=''U'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'12 = Chpgv, with ITYPE=3 and UPLO=''U'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'13 = Chegv, with ITYPE=3 and UPLO=''L'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ',/,"
                          "'14 = Chpgv, with ITYPE=3 and UPLO=''L'':',"
                          "'  | B A Z - Z D | / ( |A| |Z| n ulp )     ')");
            //
        }
        //
    } else if (Mlsamen(2, c2, "BD")) {
        //
        if (sord) {
            //
            //           Real Singular Value Decomposition:
            //
            write(iounit, "(/,1x,a3,' -- Real Singular Value Decomposition')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9973);
            //
            //           Tests performed
            //
            write(iounit, format_9972), "orthogonal";
            write(iounit, format_9971);
        } else {
            //
            //           Complex Singular Value Decomposition:
            //
            write(iounit, "(/,1x,a3,' -- Complex Singular Value Decomposition')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9973);
            //
            //           Tests performed
            //
            write(iounit, format_9972), "unitary   ";
            write(iounit, format_9971);
        }
        //
    } else if (Mlsamen(2, c2, "BB")) {
        //
        if (sord) {
            //
            //           Real General Band reduction to bidiagonal form:
            //
            write(iounit, "(/,1x,a3,' -- Real Band reduc. to bidiagonal form')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9970);
            //
            //           Tests performed
            //
            write(iounit, format_9969), "orthogonal";
        } else {
            //
            //           Complex Band reduction to bidiagonal form:
            //
            write(iounit, "(/,1x,a3,' -- Complex Band reduc. to bidiagonal form')"), path;
            //
            //           Matrix types
            //
            write(iounit, format_9970);
            //
            //           Tests performed
            //
            write(iounit, format_9969), "unitary   ";
        }
        //
    } else {
        //
        write(iounit, format_9999), path;
        return;
    }
    //
    //     Symmetric/Hermitian eigenproblem
    //
    //     Symmetric/Hermitian Generalized eigenproblem
    //
    //     Singular Value Decomposition
    //
    //     Band reduction to bidiagonal form
    //
    //     End of Rlahd2
    //
}
