--- a/TESTING/LIN/zerrtsqr.f	2024-12-03 20:39:11.000000000 +0900
+++ b/TESTING/LIN/zerrtsqr.f	2026-03-26 17:49:23.713365471 +0900
@@ -73,8 +73,8 @@
       INTEGER            I, INFO, J, MB, NB
 *     ..
 *     .. Local Arrays ..
-      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ),
-     $                   C( NMAX, NMAX ), TAU(NMAX)
+      COMPLEX*16         A( NMAX, NMAX ), W( NMAX ),
+     $                   C( NMAX, NMAX ), TAU( 5 )
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           ALAESM, CHKXER, ZGEQR,
@@ -103,7 +103,6 @@
          DO I = 1, NMAX
             A( I, J ) = 1.D0 / DBLE( I+J )
             C( I, J ) = 1.D0 / DBLE( I+J )
-            T( I, J ) = 1.D0 / DBLE( I+J )
          END DO
          W( J ) = 0.D0
       END DO
@@ -161,8 +160,9 @@
 *
 *     ZGEMQR
 *
-      TAU(1)=1
-      TAU(2)=1
+      DO I = 1, 5
+         TAU( I ) = 1
+      END DO
       SRNAMT = 'ZGEMQR'
       NB=1
       INFOT = 1
@@ -251,8 +251,9 @@
 *
 *     ZGEMLQ
 *
-      TAU(1)=1
-      TAU(2)=1
+      DO I = 1, 5
+         TAU( I ) = 1
+      END DO
       SRNAMT = 'ZGEMLQ'
       NB=1
       INFOT = 1
--- a/TESTING/LIN/derrtsqr.f	2024-12-03 20:39:11.000000000 +0900
+++ b/TESTING/LIN/derrtsqr.f	2026-03-23 12:53:49.285978304 +0900
@@ -73,8 +73,8 @@
       INTEGER            I, INFO, J, MB, NB
 *     ..
 *     .. Local Arrays ..
-      DOUBLE PRECISION   A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ),
-     $                   C( NMAX, NMAX ), TAU(NMAX*2)
+      DOUBLE PRECISION   A( NMAX, NMAX ), W( NMAX ),
+     $                   C( NMAX, NMAX ), TAU( 5 )
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           ALAESM, CHKXER, DGEQR,
@@ -103,7 +103,6 @@
          DO I = 1, NMAX
             A( I, J ) = 1.D0 / DBLE( I+J )
             C( I, J ) = 1.D0 / DBLE( I+J )
-            T( I, J ) = 1.D0 / DBLE( I+J )
          END DO
          W( J ) = 0.D0
       END DO
@@ -161,10 +160,9 @@
 *
 *     DGEMQR
 *
-      TAU(1)=1
-      TAU(2)=1
-      TAU(3)=1
-      TAU(4)=1
+      DO I = 1, 5
+         TAU( I ) = 1
+      END DO
       SRNAMT = 'DGEMQR'
       NB=1
       INFOT = 1
@@ -253,8 +251,9 @@
 *
 *     DGEMLQ
 *
-      TAU(1)=1
-      TAU(2)=1
+      DO I = 1, 5
+         TAU( I ) = 1
+      END DO
       SRNAMT = 'DGEMLQ'
       NB=1
       INFOT = 1
