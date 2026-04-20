--- a/TESTING/EIG/ddrvst2stg.f_	2026-01-23 12:05:47.429719592 +0900
+++ b/TESTING/EIG/ddrvst2stg.f	2026-01-23 12:06:03.643106479 +0900
@@ -776,10 +776,10 @@
             IF( JTYPE.LE.7 ) THEN
                NTEST = 1
                DO 120 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   120          CONTINUE
                DO 130 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   130          CONTINUE
                SRNAMT = 'DSTEV'
                CALL DSTEV( 'V', N, D1, D2, Z, LDU, WORK, IINFO )
@@ -800,17 +800,17 @@
 *              Do tests 1 and 2.
 *
                DO 140 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   140          CONTINUE
                DO 150 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   150          CONTINUE
                CALL DSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK,
      $                      RESULT( 1 ) )
 *
                NTEST = 3
                DO 160 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   160          CONTINUE
                SRNAMT = 'DSTEV'
                CALL DSTEV( 'N', N, D3, D4, Z, LDU, WORK, IINFO )
@@ -842,10 +842,10 @@
                NTEST = 4
                DO 190 I = 1, N
                   EVEIGS( I ) = D3( I )
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   190          CONTINUE
                DO 200 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   200          CONTINUE
                SRNAMT = 'DSTEVX'
                CALL DSTEVX( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL,
@@ -873,17 +873,17 @@
 *              Do tests 4 and 5.
 *
                DO 210 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   210          CONTINUE
                DO 220 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   220          CONTINUE
                CALL DSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK,
      $                      RESULT( 4 ) )
 *
                NTEST = 6
                DO 230 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   230          CONTINUE
                SRNAMT = 'DSTEVX'
                CALL DSTEVX( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL,
@@ -917,10 +917,10 @@
 *
                NTEST = 7
                DO 260 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   260          CONTINUE
                DO 270 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   270          CONTINUE
                SRNAMT = 'DSTEVR'
                CALL DSTEVR( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL,
@@ -947,17 +947,17 @@
 *              Do tests 7 and 8.
 *
                DO 280 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   280          CONTINUE
                DO 290 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   290          CONTINUE
                CALL DSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK,
      $                      RESULT( 7 ) )
 *
                NTEST = 9
                DO 300 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   300          CONTINUE
                SRNAMT = 'DSTEVR'
                CALL DSTEVR( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL,
@@ -992,10 +992,10 @@
 *
                NTEST = 10
                DO 330 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   330          CONTINUE
                DO 340 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   340          CONTINUE
                SRNAMT = 'DSTEVX'
                CALL DSTEVX( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL,
@@ -1018,10 +1018,10 @@
 *              Do tests 10 and 11.
 *
                DO 350 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   350          CONTINUE
                DO 360 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   360          CONTINUE
                CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
      $                      MAX( 1, M2 ), RESULT( 10 ) )
@@ -1029,7 +1029,7 @@
 *
                NTEST = 12
                DO 370 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   370          CONTINUE
                SRNAMT = 'DSTEVX'
                CALL DSTEVX( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL,
@@ -1079,10 +1079,10 @@
                END IF
 *
                DO 390 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   390          CONTINUE
                DO 400 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   400          CONTINUE
                SRNAMT = 'DSTEVX'
                CALL DSTEVX( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL,
@@ -1112,17 +1112,17 @@
 *              Do tests 13 and 14.
 *
                DO 410 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   410          CONTINUE
                DO 420 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   420          CONTINUE
                CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
      $                      MAX( 1, M2 ), RESULT( 13 ) )
 *
                NTEST = 15
                DO 430 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   430          CONTINUE
                SRNAMT = 'DSTEVX'
                CALL DSTEVX( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL,
@@ -1150,10 +1150,10 @@
 *
                NTEST = 16
                DO 450 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   450          CONTINUE
                DO 460 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   460          CONTINUE
                SRNAMT = 'DSTEVD'
                CALL DSTEVD( 'V', N, D1, D2, Z, LDU, WORK, LWEDC, IWORK,
@@ -1175,17 +1175,17 @@
 *              Do tests 16 and 17.
 *
                DO 470 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   470          CONTINUE
                DO 480 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   480          CONTINUE
                CALL DSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK,
      $                      RESULT( 16 ) )
 *
                NTEST = 18
                DO 490 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   490          CONTINUE
                SRNAMT = 'DSTEVD'
                CALL DSTEVD( 'N', N, D3, D4, Z, LDU, WORK, LWEDC, IWORK,
@@ -1218,10 +1218,10 @@
 *
                NTEST = 19
                DO 520 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   520          CONTINUE
                DO 530 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   530          CONTINUE
                SRNAMT = 'DSTEVR'
                CALL DSTEVR( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL,
@@ -1244,10 +1244,10 @@
 *              DO tests 19 and 20.
 *
                DO 540 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   540          CONTINUE
                DO 550 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   550          CONTINUE
                CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
      $                      MAX( 1, M2 ), RESULT( 19 ) )
@@ -1255,7 +1255,7 @@
 *
                NTEST = 21
                DO 560 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   560          CONTINUE
                SRNAMT = 'DSTEVR'
                CALL DSTEVR( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL,
@@ -1305,10 +1305,10 @@
                END IF
 *
                DO 580 I = 1, N
-                  D1( I ) = DBLE( A( I, I ) )
+                  D1( I ) = A( I, I )
   580          CONTINUE
                DO 590 I = 1, N - 1
-                  D2( I ) = DBLE( A( I+1, I ) )
+                  D2( I ) = A( I+1, I )
   590          CONTINUE
                SRNAMT = 'DSTEVR'
                CALL DSTEVR( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL,
@@ -1338,17 +1338,17 @@
 *              Do tests 22 and 23.
 *
                DO 600 I = 1, N
-                  D3( I ) = DBLE( A( I, I ) )
+                  D3( I ) = A( I, I )
   600          CONTINUE
                DO 610 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   610          CONTINUE
                CALL DSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
      $                      MAX( 1, M2 ), RESULT( 22 ) )
 *
                NTEST = 24
                DO 620 I = 1, N - 1
-                  D4( I ) = DBLE( A( I+1, I ) )
+                  D4( I ) = A( I+1, I )
   620          CONTINUE
                SRNAMT = 'DSTEVR'
                CALL DSTEVR( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL,
