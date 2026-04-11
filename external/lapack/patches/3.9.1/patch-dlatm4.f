--- a/TESTING/EIG/dlatm4.f~	2021-03-26 03:25:15.000000000 +0900
+++ b/TESTING/EIG/dlatm4.f	2026-01-23 11:27:57.289505650 +0900
@@ -342,10 +342,10 @@
 *        Scale by AMAGN
 *
          DO 230 JD = KBEG, KEND
-            A( JD, JD ) = AMAGN*DBLE( A( JD, JD ) )
+            A( JD, JD ) = AMAGN* A( JD, JD )
   230    CONTINUE
          DO 240 JD = ISDB, ISDE
-            A( JD+1, JD ) = AMAGN*DBLE( A( JD+1, JD ) )
+            A( JD+1, JD ) = AMAGN* A( JD+1, JD )
   240    CONTINUE
 *
 *        If ISIGN = 1 or 2, assign random signs to diagonal and
@@ -353,13 +353,13 @@
 *
          IF( ISIGN.GT.0 ) THEN
             DO 250 JD = KBEG, KEND
-               IF( DBLE( A( JD, JD ) ).NE.ZERO ) THEN
+               IF( A( JD, JD ).NE.ZERO ) THEN
                   IF( DLARAN( ISEED ).GT.HALF )
      $               A( JD, JD ) = -A( JD, JD )
                END IF
   250       CONTINUE
             DO 260 JD = ISDB, ISDE
-               IF( DBLE( A( JD+1, JD ) ).NE.ZERO ) THEN
+               IF( A( JD+1, JD ).NE.ZERO ) THEN
                   IF( DLARAN( ISEED ).GT.HALF )
      $               A( JD+1, JD ) = -A( JD+1, JD )
                END IF
