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

// Derived from LAPACK routine ALAHD.
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

void Alahd(INTEGER const iounit, fem::str_cref path) {
    common cmn;
    common_write write(cmn);
    //
    // First line of header
    //
    static const char *format_9999 = "(/,1x,a3,':  General dense matrices')";
    static const char *format_9998 = "(/,1x,a3,':  General band matrices')";
    static const char *format_9997 = "(/,1x,a3,':  General tridiagonal')";
    static const char *format_9996 = "(/,1x,a3,':  ',a9,' positive definite matrices')";
    static const char *format_9995 = "(/,1x,a3,':  ',a9,' positive definite packed matrices')";
    static const char *format_9994 = "(/,1x,a3,':  ',a9,' positive definite band matrices')";
    static const char *format_9993 = "(/,1x,a3,':  ',a9,' positive definite tridiagonal')";
    static const char *format_9992 = "(/,1x,a3,':  ',a9,' indefinite matrices',"
                                     "', partial (Bunch-Kaufman) pivoting')";
    static const char *format_9991 = "(/,1x,a3,':  ',a9,' indefinite packed matrices',"
                                     "', partial (Bunch-Kaufman) pivoting')";
    static const char *format_9892 = "(/,1x,a3,':  ',a9,' indefinite matrices',"
                                     "', \"rook\" (bounded Bunch-Kaufman) pivoting')";
    static const char *format_9990 = "(/,1x,a3,':  Triangular matrices')";
    static const char *format_9989 = "(/,1x,a3,':  Triangular packed matrices')";
    static const char *format_9988 = "(/,1x,a3,':  Triangular band matrices')";
    static const char *format_9987 = "(/,1x,a3,':  ',a2,' factorization of general matrices')";
    static const char *format_9985 = "(/,1x,a3,':  RQ factorization of trapezoidal matrix')";
    static const char *format_9984 = "(/,1x,a3,':  Least squares driver routines')";
    static const char *format_9983 = "(/,1x,a3,':  LU factorization variants')";
    static const char *format_9982 = "(/,1x,a3,':  Cholesky factorization variants')";
    static const char *format_9981 = "(/,1x,a3,':  QR factorization variants')";
    static const char *format_9980 = "(/,1x,a3,':  No header available')";
    static const char *format_8000 = "(/,1x,a3,':  QRT factorization for general matrices')";
    static const char *format_8001 = "(/,1x,a3,':  QRT factorization for ','triangular-pentagonal matrices')";
    static const char *format_8002 = "(/,1x,a3,':  LQT factorization for general matrices')";
    static const char *format_8003 = "(/,1x,a3,':  LQT factorization for ','triangular-pentagonal matrices')";
    static const char *format_8004 = "(/,1x,a3,':  TS factorization for ','tall-skinny or short-wide matrices')";
    static const char *format_8005 = "(/,1x,a3,':  Householder reconstruction from TSQR',"
                                     "' factorization output ',/,' for tall-skinny matrices.')";
    static const char *format_8006 = "(/,1x,a3,':  truncated QR factorization',' with column pivoting')";
    //
    // GE matrix types
    //
    static const char *format_9979 = "(4x,'1. Diagonal',24x,'7. Last n/2 columns zero',/,4x,"
                                     "'2. Upper triangular',16x,'8. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. Lower triangular',16x,'9. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Random, CNDNUM = 2',13x,'10. Scaled near underflow',/,4x,"
                                     "'5. First column zero',14x,'11. Scaled near overflow',/,4x,"
                                     "'6. Last column zero')";
    //
    // GB matrix types
    //
    static const char *format_9978 = "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. First column zero',15x,'6. Random, CNDNUM = .01/EPS',/,4x,"
                                     "'3. Last column zero',16x,'7. Scaled near underflow',/,4x,"
                                     "'4. Last n/2 columns zero',11x,'8. Scaled near overflow')";
    //
    // GT matrix types
    //
    static const char *format_9977 = "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                                     "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. First column zero',/,4x,"
                                     "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last column zero',/,4x,"
                                     "'4. Random, CNDNUM = 0.1/EPS',7x,'10. Last n/2 columns zero',/,4x,"
                                     "'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                                     "'6. Scaled near overflow',11x,'12. Scaled near overflow')";
    //
    // PT matrix types
    //
    static const char *format_9976 = "(' Matrix types (1-6 have specified condition numbers):',/,4x,"
                                     "'1. Diagonal',24x,'7. Random, unspecified CNDNUM',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. First row and column zero',/,4x,"
                                     "'3. Random, CNDNUM = sqrt(0.1/EPS)',2x,'9. Last row and column zero',/,"
                                     "4x,'4. Random, CNDNUM = 0.1/EPS',7x,'10. Middle row and column zero',/,"
                                     "4x,'5. Scaled near underflow',10x,'11. Scaled near underflow',/,4x,"
                                     "'6. Scaled near overflow',11x,'12. Scaled near overflow')";
    //
    // PO, PP matrix types
    //
    static const char *format_9975 = "(4x,'1. Diagonal',24x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = 0.1/EPS',/,3x,"
                                     "'*3. First row and column zero',7x,'8. Scaled near underflow',/,3x,"
                                     "'*4. Last row and column zero',8x,'9. Scaled near overflow',/,3x,"
                                     "'*5. Middle row and column zero',/,3x,'(* - tests error exits from ',a3,"
                                     "'TRF, no test ratios are computed)')";
    //
    // CH matrix types
    //
    static const char *format_9974 = "(4x,'1. Diagonal',24x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = 0.1/EPS',/,3x,"
                                     "'*3. First row and column zero',7x,'8. Scaled near underflow',/,3x,"
                                     "'*4. Last row and column zero',8x,'9. Scaled near overflow',/,3x,"
                                     "'*5. Middle row and column zero',/,3x,"
                                     "'(* - tests error exits, no test ratios are computed)')";
    //
    // PS matrix types
    //
    static const char *format_8973 = "(4x,'1. Diagonal',/,4x,'2. Random, CNDNUM = 2',14x,/,3x,"
                                     "'*3. Nonzero eigenvalues of: D(1:RANK-1)=1 and ','D(RANK) = 1.0/',a4,/,"
                                     "3x,'*4. Nonzero eigenvalues of: D(1)=1 and ',' D(2:RANK) = 1.0/',a4,/,3x,"
                                     "'*5. Nonzero eigenvalues of: D(I) = ',a4,'**(-(I-1)/(RANK-1)) ',"
                                     "' I=1:RANK',/,4x,'6. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'7. Random, CNDNUM = 0.1/EPS',/,4x,'8. Scaled near underflow',/,4x,"
                                     "'9. Scaled near overflow',/,3x,'(* - Semi-definite tests )')";
    static const char *format_8972 = "(3x,'RANK minus computed rank, returned by ',a,'PSTRF')";
    //
    // PB matrix types
    //
    static const char *format_9973 = "(4x,'1. Random, CNDNUM = 2',14x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,3x,"
                                     "'*2. First row and column zero',7x,'6. Random, CNDNUM = 0.1/EPS',/,3x,"
                                     "'*3. Last row and column zero',8x,'7. Scaled near underflow',/,3x,"
                                     "'*4. Middle row and column zero',6x,'8. Scaled near overflow',/,3x,"
                                     "'(* - tests error exits from ',a3,'TRF, no test ratios are computed)')";
    //
    // SSY, SSR, SSP, CHE, CHR, CHP matrix types
    //
    static const char *format_9972 = "(4x,'1. Diagonal',24x,'6. Last n/2 rows and columns zero',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'3. First row and column zero',7x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'4. Last row and column zero',8x,'9. Scaled near underflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'10. Scaled near overflow')";
    //
    // CSY, CSR, CSP matrix types
    //
    static const char *format_9971 = "(4x,'1. Diagonal',24x,'7. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Random, CNDNUM = 2',14x,'8. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. First row and column zero',7x,'9. Scaled near underflow',/,4x,"
                                     "'4. Last row and column zero',7x,'10. Scaled near overflow',/,4x,"
                                     "'5. Middle row and column zero',5x,'11. Block diagonal matrix',/,4x,"
                                     "'6. Last n/2 rows and columns zero')";
    //
    // QR matrix types
    //
    static const char *format_9970 = "(4x,'1. Diagonal',24x,'5. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'2. Upper triangular',16x,'6. Random, CNDNUM = 0.1/EPS',/,4x,"
                                     "'3. Lower triangular',16x,'7. Scaled near underflow',/,4x,"
                                     "'4. Random, CNDNUM = 2',14x,'8. Scaled near overflow')";
    //
    // QP matrix types
    //
    static const char *format_9969 = "(' Matrix types (2-6 have condition 1/EPS):',/,4x,'1. Zero matrix',21x,"
                                     "'4. First n/2 columns fixed',/,4x,'2. One small eigenvalue',12x,"
                                     "'5. Last n/2 columns fixed',/,4x,'3. Geometric distribution',10x,"
                                     "'6. Every second column fixed')";
    //
    // QK matrix types
    //
    static const char *format_9871 = "(4x,' 1. Zero matrix',/,4x,' 2. Random, Diagonal, CNDNUM = 2',/,4x,"
                                     "' 3. Random, Upper triangular, CNDNUM = 2',/,4x,"
                                     "' 4. Random, Lower triangular, CNDNUM = 2',/,4x,"
                                     "' 5. Random, First column is zero, CNDNUM = 2',/,4x,"
                                     "' 6. Random, Last MINMN column is zero, CNDNUM = 2',/,4x,"
                                     "' 7. Random, Last N column is zero, CNDNUM = 2',/,4x,"
                                     "' 8. Random, Middle column in MINMN is zero,',' CNDNUM = 2',/,4x,"
                                     "' 9. Random, First half of MINMN columns are zero,',' CNDNUM = 2',/,4x,"
                                     "'10. Random, Last columns are zero starting from',"
                                     "' MINMN/2+1, CNDNUM = 2',/,4x,"
                                     "'11. Random, Half MINMN columns in the middle are',"
                                     "' zero starting from MINMN/2-(MINMN/2)/2+1,',' CNDNUM = 2',/,4x,"
                                     "'12. Random, Odd columns are ZERO, CNDNUM = 2',/,4x,"
                                     "'13. Random, Even columns are ZERO, CNDNUM = 2',/,4x,"
                                     "'14. Random, CNDNUM = 2',/,4x,'15. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                     "'16. Random, CNDNUM = 0.1/EPS',/,4x,'17. Random, CNDNUM = 0.1/EPS,',"
                                     "' one small singular value S(N)=1/CNDNUM',/,4x,"
                                     "'18. Random, CNDNUM = 2, scaled near underflow,',"
                                     "' NORM = SMALL = SAFMIN',/,4x,"
                                     "'19. Random, CNDNUM = 2, scaled near overflow,',"
                                     "' NORM = LARGE = 1.0/( 0.25 * ( SAFMIN / EPS ) )')";
    //
    // TZ matrix types
    //
    static const char *format_9968 = "(' Matrix types (2-3 have condition 1/EPS):',/,4x,'1. Zero matrix',/,4x,"
                                     "'2. One small eigenvalue',/,4x,'3. Geometric distribution')";
    //
    // LS matrix types
    //
    static const char *format_9967 = "(' Matrix types (1-3: full rank, 4-6: rank deficient):',/,4x,"
                                     "'1 and 4. Normal scaling',/,4x,'2 and 5. Scaled near overflow',/,4x,"
                                     "'3 and 6. Scaled near underflow')";
    //
    // TR, TP matrix types
    //
    static const char *format_9966 = "(' Matrix types for ',a3,' routines:',/,4x,'1. Diagonal',24x,"
                                     "'6. Scaled near overflow',/,4x,'2. Random, CNDNUM = 2',14x,'7. Identity',"
                                     "/,4x,'3. Random, CNDNUM = sqrt(0.1/EPS)  ',"
                                     "'8. Unit triangular, CNDNUM = 2',/,4x,'4. Random, CNDNUM = 0.1/EPS',8x,"
                                     "'9. Unit, CNDNUM = sqrt(0.1/EPS)',/,4x,'5. Scaled near underflow',10x,"
                                     "'10. Unit, CNDNUM = 0.1/EPS')";
    static const char *format_9965 = "(' Special types for testing ',a,':',/,3x,"
                                     "'11. Matrix elements are O(1), large right hand side',/,3x,"
                                     "'12. First diagonal causes overflow,',' offdiagonal column norms < 1',/,"
                                     "3x,'13. First diagonal causes overflow,',' offdiagonal column norms > 1',"
                                     "/,3x,'14. Growth factor underflows, solution does not overflow',/,3x,"
                                     "'15. Small diagonal causes gradual overflow',/,3x,"
                                     "'16. One zero diagonal element',/,3x,"
                                     "'17. Large offdiagonals cause overflow when adding a column',/,3x,"
                                     "'18. Unit triangular with large right hand side')";
    //
    // TB matrix types
    //
    static const char *format_9964 = "(' Matrix types for ',a3,' routines:',/,4x,'1. Random, CNDNUM = 2',14x,"
                                     "'6. Identity',/,4x,'2. Random, CNDNUM = sqrt(0.1/EPS)  ',"
                                     "'7. Unit triangular, CNDNUM = 2',/,4x,'3. Random, CNDNUM = 0.1/EPS',8x,"
                                     "'8. Unit, CNDNUM = sqrt(0.1/EPS)',/,4x,'4. Scaled near underflow',11x,"
                                     "'9. Unit, CNDNUM = 0.1/EPS',/,4x,'5. Scaled near overflow')";
    static const char *format_9963 = "(' Special types for testing ',a,':',/,3x,"
                                     "'10. Matrix elements are O(1), large right hand side',/,3x,"
                                     "'11. First diagonal causes overflow,',' offdiagonal column norms < 1',/,"
                                     "3x,'12. First diagonal causes overflow,',' offdiagonal column norms > 1',"
                                     "/,3x,'13. Growth factor underflows, solution does not overflow',/,3x,"
                                     "'14. Small diagonal causes gradual overflow',/,3x,"
                                     "'15. One zero diagonal element',/,3x,"
                                     "'16. Large offdiagonals cause overflow when adding a column',/,3x,"
                                     "'17. Unit triangular with large right hand side')";
    //
    // Test ratios
    //
    static const char *format_9962 = "(3x,i2,': norm( L * U - A )  / ( N * norm(A) * EPS )')";
    static const char *format_9961 = "(3x,i2,': norm( I - A*AINV ) / ','( N * norm(A) * norm(AINV) * EPS )')";
    static const char *format_9960 = "(3x,i2,': norm( B - A * X )  / ','( norm(A) * norm(X) * EPS )')";
    static const char *format_6660 = "(3x,i2,': diagonal is not non-negative')";
    static const char *format_9959 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * CNDNUM * EPS )')";
    static const char *format_9958 = "(3x,i2,': norm( X - XACT )   / ',"
                                     "'( norm(XACT) * CNDNUM * EPS ), refined')";
    static const char *format_9957 = "(3x,i2,': norm( X - XACT )   / ','( norm(XACT) * (error bound) )')";
    static const char *format_9956 = "(3x,i2,': (backward error)   / EPS')";
    static const char *format_9955 = "(3x,i2,': RCOND * CNDNUM - 1.0')";
    static const char *format_9954 = "(3x,i2,': norm( U'' * U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L * L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_8950 = "(3x,'norm( P * U'' * U * P'' - A ) / ( N * norm(A) * EPS )',', or',/,3x,"
                                     "'norm( P * L * L'' * P'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9953 = "(3x,i2,': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9952 = "(3x,i2,': norm( U''*D*U - A ) / ( N * norm(A) * EPS )',', or',/,7x,"
                                     "'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')";
    static const char *format_9951 = "(' Test ratio for ',a,':',/,3x,i2,"
                                     "': norm( s*b - A*x )  / ( norm(A) * norm(x) * EPS )')";
    static const char *format_9950 = "(3x,i2,': norm( R - Q'' * A ) / ( M * norm(A) * EPS )')";
    static const char *format_6950 = "(3x,i2,': norm( R - Q'' * A ) / ( M * norm(A) * EPS )       [RFPG]')";
    static const char *format_9949 = "(3x,i2,': norm( L - A * Q'' ) / ( N * norm(A) * EPS )')";
    static const char *format_9948 = "(3x,i2,': norm( L - Q'' * A ) / ( M * norm(A) * EPS )')";
    static const char *format_9947 = "(3x,i2,': norm( R - A * Q'' ) / ( N * norm(A) * EPS )')";
    static const char *format_9946 = "(3x,i2,': norm( I - Q''*Q )   / ( M * EPS )')";
    static const char *format_9945 = "(3x,i2,': norm( I - Q*Q'' )   / ( N * EPS )')";
    static const char *format_9944 = "(3x,i2,': norm( Q*C - Q*C )  / ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9943 = "(3x,i2,': norm( C*Q - C*Q )  / ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9942 = "(3x,i2,': norm( Q''*C - Q''*C )/ ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9941 = "(3x,i2,': norm( C*Q'' - C*Q'' )/ ','( ',a1,' * norm(C) * EPS )')";
    static const char *format_9940 = "(3x,i2,': norm(svd(A) - svd(R)) / ','( M * norm(svd(R)) * EPS )')";
    static const char *format_9939 = "(3x,i2,': norm( A*P - Q*R ) / ( M * norm(A) * EPS )')";
    static const char *format_9938 = "(3x,i2,': norm( I - Q''*Q ) / ( M * EPS )')";
    static const char *format_9937 = "(3x,i2,': norm( A - R*Q )       / ( M * norm(A) * EPS )')";
    static const char *format_9935 = "(3x,i2,': norm( B - A * X )   / ',"
                                     "'( max(M,N) * norm(A) * norm(X) * EPS )')";
    static const char *format_9934 = "(3x,i2,': norm( (A*X-B)'' *A ) / ',"
                                     "'( max(M,N,NRHS) * norm(A) * norm(B) * EPS )')";
    static const char *format_9933 = "(3x,i2,': norm(svd(A)-svd(R)) / ','( min(M,N) * norm(svd(R)) * EPS )')";
    static const char *format_9932 = "(3x,i2,': Check if X is in the row space of A or A''')";
    static const char *format_9931 = "(3x,i2,': norm( (A*X-B)'' *A ) / ',"
                                     "'( max(M,N,NRHS) * norm(A) * norm(B) * EPS )',/,7x,"
                                     "'if TRANS=''N'' and M.GE.N or TRANS=''T'' and M.LT.N, ','otherwise',/,7x,"
                                     "'check if X is in the row space of A or A'' ','(overdetermined case)')";
    static const char *format_9929 = "(' Test ratios (1-3: ',a1,'TZRZF):')";
    static const char *format_9919 = "(3x,' 3-4: same as 1-2',3x,' 5-6: same as 1-2')";
    static const char *format_9920 = "(3x,' 11-14: same as 7-10',3x,' 15-18: same as 7-10')";
    static const char *format_9921 = "(' Test ratios:',/,'    (1-2: ',a1,'GELS, 3-4: ',a1,'GELST, 5-6: ',a1,"
                                     "'GETSLS, 7-10: ',a1,'GELSY, 11-14: ',a1,'GETSS, 15-18: ',a1,'GELSD)')";
    static const char *format_9928 = "(7x,'where ALPHA = ( 1 + SQRT( 17 ) ) / 8')";
    static const char *format_9927 = "(3x,i2,': ABS( Largest element in L )',/,12x,"
                                     "' - ( 1 / ( 1 - ALPHA ) ) + THRESH')";
    static const char *format_9926 = "(3x,i2,': Largest 2-Norm of 2-by-2 pivots',/,12x,"
                                     "' - ( ( 1 + ALPHA ) / ( 1 - ALPHA ) ) + THRESH')";
    static const char *format_8011 = "(3x,i2,': norm( R - Q''*A ) / ( M * norm(A) * EPS )')";
    static const char *format_8012 = "(3x,i2,': norm( I - Q''*Q ) / ( M * EPS )')";
    static const char *format_8013 = "(3x,i2,': norm( Q*C - Q*C ) / ( M * norm(C) * EPS )')";
    static const char *format_8014 = "(3x,i2,': norm( Q''*C - Q''*C ) / ( M * norm(C) * EPS )')";
    static const char *format_8015 = "(3x,i2,': norm( C*Q - C*Q ) / ( M * norm(C) * EPS )')";
    static const char *format_8016 = "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( M * norm(C) * EPS )')";
    static const char *format_8017 = "(3x,i2,': norm( R - Q''*A ) / ( (M+N) * norm(A) * EPS )')";
    static const char *format_8018 = "(3x,i2,': norm( I - Q''*Q ) / ( (M+N) * EPS )')";
    static const char *format_8019 = "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8020 = "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8021 = "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8022 = "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8023 = "(3x,i2,': norm( L - A*Q'' ) / ( (M+N) * norm(A) * EPS )')";
    static const char *format_8024 = "(3x,i2,': norm( I - Q*Q'' ) / ( (M+N) * EPS )')";
    static const char *format_8025 = "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8026 = "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8027 = "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8028 = "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8029 = "(3x,i2,': norm( L - A*Q'' ) / ( (M+N) * norm(A) * EPS )')";
    static const char *format_8030 = "(3x,i2,': norm( I - Q*Q'' ) / ( (M+N) * EPS )')";
    static const char *format_8031 = "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8032 = "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8033 = "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8034 = "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8035 = "(3x,i2,': norm( R - Q''*A ) / ( (M+N) * norm(A) * EPS )')";
    static const char *format_8036 = "(3x,i2,': norm( I - Q''*Q ) / ( (M+N) * EPS )')";
    static const char *format_8037 = "(3x,i2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8038 = "(3x,i2,': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8039 = "(3x,i2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')";
    static const char *format_8040 = "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')";
    //
    static const char *format_8050 = "(3x,i2,': norm( R - Q''*A ) / ( M * norm(A) * EPS )')";
    static const char *format_8051 = "(3x,i2,': norm( I - Q''*Q ) / ( M * EPS )')";
    static const char *format_8052 = "(3x,i2,': norm( Q*C - Q*C ) / ( M * norm(C) * EPS )')";
    static const char *format_8053 = "(3x,i2,': norm( Q''*C - Q''*C ) / ( M * norm(C) * EPS )')";
    static const char *format_8054 = "(3x,i2,': norm( C*Q - C*Q ) / ( M * norm(C) * EPS )')";
    static const char *format_8055 = "(3x,i2,': norm( C*Q'' - C*Q'' ) / ( M * norm(C) * EPS )')";
    //
    static const char *format_8060 = "(3x,i2,': 2-norm(svd(A) - svd(R)) / ',"
                                     "'( max(M,N) * 2-norm(svd(R)) * EPS )')";
    static const char *format_8061 = "(3x,i2,': 1-norm( A*P - Q*R ) / ( max(M,N) * 1-norm(A)',' * EPS )')";
    static const char *format_8062 = "(3x,i2,': 1-norm( I - Q''*Q ) / ( M * EPS )')";
    static const char *format_8063 = "(3x,i2,': Returns 1.0D+100, if abs(R(K+1,K+1))',"
                                     "' > abs(R(K,K)), where K=1:KFACT-1')";
    static const char *format_8064 = "(3x,i2,': 1-norm(Q**T * B - Q**T * B ) / ( M * EPS )')";
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
    fem::str<4> eigcnm;
    fem::str<32> subnam;
    if (Mlsamen(2, p2.elems, "GE")) {
        //
        // GE: General dense
        //
        write(iounit, format_9999), path;
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
    } else if (Mlsamen(2, p2.elems, "GB")) {
        //
        // GB: General band
        //
        write(iounit, format_9998), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9978);
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
    } else if (Mlsamen(2, p2.elems, "GT")) {
        //
        // GT: General tridiagonal
        //
        write(iounit, format_9997), path;
        write(iounit, format_9977);
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
    } else if (Mlsamen(2, p2.elems, "PO") || Mlsamen(2, p2.elems, "PP")) {
        //
        // PO: Positive definite full
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
        write(iounit, format_9975), path;
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
    } else if (Mlsamen(2, p2.elems, "PS")) {
        //
        // PS: Positive semi-definite full
        //
        if (sord) {
            sym = "Symmetric";
        } else {
            sym = "Hermitian";
        }
        if (Mlsame(c1.elems, "S") || Mlsame(c1.elems, "C")) {
            eigcnm = "1E04";
        } else {
            eigcnm = "1D12";
        }
        write(iounit, format_9995), path, sym;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_8973), eigcnm, eigcnm, eigcnm;
        write(iounit, "(' Difference:')");
        write(iounit, format_8972), c1;
        write(iounit, "(' Test ratio:')");
        write(iounit, format_8950);
        write(iounit, "(' Messages:')");
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
        write(iounit, format_9973), path;
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
    } else if (Mlsamen(2, p2.elems, "PT")) {
        //
        // PT: Positive definite tridiagonal
        //
        if (sord) {
            write(iounit, format_9993), path, "Symmetric";
        } else {
            write(iounit, format_9993), path, "Hermitian";
        }
        write(iounit, format_9976);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9952), 1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "SY")) {
        //
        // SY: Symmetric indefinite full,
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3.elems, "Y")) {
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
    } else if (Mlsamen(2, p2.elems, "SP")) {
        //
        // SP: Symmetric indefinite packed,
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3.elems, "Y")) {
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
    } else if (Mlsamen(2, p2.elems, "HA")) {
        //
        // HA: Hermitian,
        // with Assen Algorithm
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
    } else if (Mlsamen(2, p2.elems, "HE")) {
        //
        // HE: Hermitian indefinite full,
        // with partial (Bunch-Kaufman) pivoting algorithm
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
    } else if (Mlsamen(2, p2.elems, "HR") || Mlsamen(2, p2.elems, "HR")) {
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
    } else if (Mlsamen(2, p2.elems, "HP")) {
        //
        // HP: Hermitian indefinite packed,
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        if (Mlsame(c3.elems, "E")) {
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
    } else if (Mlsamen(2, p2.elems, "TR") || Mlsamen(2, p2.elems, "TP")) {
        //
        // TR: Triangular full
        // TP: Triangular packed
        //
        if (Mlsame(c3.elems, "R")) {
            write(iounit, format_9990), path;
            subnam = path(1, 1) + "LATRS";
        } else {
            write(iounit, format_9989), path;
            subnam = path(1, 1) + "LATPS";
        }
        write(iounit, format_9966), path;
        write(iounit, format_9965), subnam(1, fem::len_trim(subnam));
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9961), 1;
        write(iounit, format_9960), 2;
        write(iounit, format_9959), 3;
        write(iounit, format_9958), 4;
        write(iounit, format_9957), 5;
        write(iounit, format_9956), 6;
        write(iounit, format_9955), 7;
        write(iounit, format_9951), subnam(1, fem::len_trim(subnam)), 8;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "TB")) {
        //
        // TB: Triangular band
        //
        write(iounit, format_9988), path;
        subnam = path(1, 1) + "LATBS";
        write(iounit, format_9964), path;
        write(iounit, format_9963), subnam(1, fem::len_trim(subnam));
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9960), 1;
        write(iounit, format_9959), 2;
        write(iounit, format_9958), 3;
        write(iounit, format_9957), 4;
        write(iounit, format_9956), 5;
        write(iounit, format_9955), 6;
        write(iounit, format_9951), subnam(1, fem::len_trim(subnam)), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "QR")) {
        //
        // QR decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "QR";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9950), 1;
        write(iounit, format_6950), 8;
        write(iounit, format_9946), 2;
        write(iounit, format_9944), 3, "M";
        write(iounit, format_9943), 4, "M";
        write(iounit, format_9942), 5, "M";
        write(iounit, format_9941), 6, "M";
        write(iounit, format_9960), 7;
        write(iounit, format_6660), 9;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "LQ")) {
        //
        // LQ decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "LQ";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9949), 1;
        write(iounit, format_9945), 2;
        write(iounit, format_9944), 3, "N";
        write(iounit, format_9943), 4, "N";
        write(iounit, format_9942), 5, "N";
        write(iounit, format_9941), 6, "N";
        write(iounit, format_9960), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "QL")) {
        //
        // QL decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "QL";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9948), 1;
        write(iounit, format_9946), 2;
        write(iounit, format_9944), 3, "M";
        write(iounit, format_9943), 4, "M";
        write(iounit, format_9942), 5, "M";
        write(iounit, format_9941), 6, "M";
        write(iounit, format_9960), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "RQ")) {
        //
        // RQ decomposition of rectangular matrices
        //
        write(iounit, format_9987), path, "RQ";
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9947), 1;
        write(iounit, format_9945), 2;
        write(iounit, format_9944), 3, "N";
        write(iounit, format_9943), 4, "N";
        write(iounit, format_9942), 5, "N";
        write(iounit, format_9941), 6, "N";
        write(iounit, format_9960), 7;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "QP")) {
        //
        // QR decomposition with column pivoting
        //
        write(iounit, format_8006), path;
        write(iounit, format_9969);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9940), 1;
        write(iounit, format_9939), 2;
        write(iounit, format_9938), 3;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "QK")) {
        //
        // truncated QR decomposition with column pivoting
        //
        write(iounit, format_8006), path;
        write(iounit, format_9871);
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8060), 1;
        write(iounit, format_8061), 2;
        write(iounit, format_8062), 3;
        write(iounit, format_8063), 4;
        write(iounit, format_8064), 5;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "TZ")) {
        //
        // TZ:  Trapezoidal
        //
        write(iounit, format_9985), path;
        write(iounit, format_9968);
        write(iounit, format_9929), c1;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_9940), 1;
        write(iounit, format_9937), 2;
        write(iounit, format_9938), 3;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "LS")) {
        //
        // LS:  Least Squares driver routines for
        // LS, LST, TSLS, LSD, LSS, LSX and LSY.
        //
        write(iounit, format_9984), path;
        write(iounit, format_9967);
        write(iounit, format_9921), c1, c1, c1, c1, c1, c1;
        write(iounit, format_9935), 1;
        write(iounit, format_9931), 2;
        write(iounit, format_9919);
        write(iounit, format_9933), 7;
        write(iounit, format_9935), 8;
        write(iounit, format_9934), 9;
        write(iounit, format_9932), 10;
        write(iounit, format_9920);
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "LU")) {
        //
        // LU factorization variants
        //
        write(iounit, format_9983), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9979);
        write(iounit, "(' Test ratio:')");
        write(iounit, format_9962), 1;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "CH")) {
        //
        // Cholesky factorization variants
        //
        write(iounit, format_9982), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9974);
        write(iounit, "(' Test ratio:')");
        write(iounit, format_9954), 1;
        write(iounit, "(' Messages:')");
        //
    } else if (Mlsamen(2, p2.elems, "QS")) {
        //
        // QR factorization variants
        //
        write(iounit, format_9981), path;
        write(iounit, "(' Matrix types:')");
        write(iounit, format_9970);
        write(iounit, "(' Test ratios:')");
        //
    } else if (Mlsamen(2, p2.elems, "QT")) {
        //
        // QRT (general matrices)
        //
        write(iounit, format_8000), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8011), 1;
        write(iounit, format_8012), 2;
        write(iounit, format_8013), 3;
        write(iounit, format_8014), 4;
        write(iounit, format_8015), 5;
        write(iounit, format_8016), 6;
        //
    } else if (Mlsamen(2, p2.elems, "QX")) {
        //
        // QRT (triangular-pentagonal)
        //
        write(iounit, format_8001), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8017), 1;
        write(iounit, format_8018), 2;
        write(iounit, format_8019), 3;
        write(iounit, format_8020), 4;
        write(iounit, format_8021), 5;
        write(iounit, format_8022), 6;
        //
    } else if (Mlsamen(2, p2.elems, "TQ")) {
        //
        // QRT (triangular-pentagonal)
        //
        write(iounit, format_8002), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8023), 1;
        write(iounit, format_8024), 2;
        write(iounit, format_8025), 3;
        write(iounit, format_8026), 4;
        write(iounit, format_8027), 5;
        write(iounit, format_8028), 6;
        //
    } else if (Mlsamen(2, p2.elems, "XQ")) {
        //
        // QRT (triangular-pentagonal)
        //
        write(iounit, format_8003), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8029), 1;
        write(iounit, format_8030), 2;
        write(iounit, format_8031), 3;
        write(iounit, format_8032), 4;
        write(iounit, format_8033), 5;
        write(iounit, format_8034), 6;
        //
    } else if (Mlsamen(2, p2.elems, "TS")) {
        //
        // TS:  QR routines for tall-skinny and short-wide matrices
        //
        write(iounit, format_8004), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8035), 1;
        write(iounit, format_8036), 2;
        write(iounit, format_8037), 3;
        write(iounit, format_8038), 4;
        write(iounit, format_8039), 5;
        write(iounit, format_8040), 6;
        //
    } else if (Mlsamen(2, p2.elems, "HH")) {
        //
        // HH:  Householder reconstruction for tall-skinny matrices
        //
        write(iounit, format_8005), path;
        write(iounit, "(' Test ratios:')");
        write(iounit, format_8050), 1;
        write(iounit, format_8051), 2;
        write(iounit, format_8052), 3;
        write(iounit, format_8053), 4;
        write(iounit, format_8054), 5;
        write(iounit, format_8055), 6;
        //
    } else {
        //
        // Print error message if no header is available.
        //
        write(iounit, format_9980), path;
    }
    //
    // End of Alahd
    //
}
