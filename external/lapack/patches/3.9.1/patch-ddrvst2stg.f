--- TESTING/EIG/ddrvst2stg.f~	2021-03-26 03:25:15.000000000 +0900
+++ TESTING/EIG/ddrvst2stg.f	2026-01-23 11:32:44.174292705 +0900
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
