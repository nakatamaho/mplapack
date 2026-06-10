diff --git a/SRC/claqp2rk.f b/SRC/claqp2rk.f
index 12ebf2d94..d27d978e9 100644
--- a/SRC/claqp2rk.f
+++ b/SRC/claqp2rk.f
@@ -358,8 +358,8 @@
       PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
 *     ..
 *     .. Local Scalars ..
-      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP, MINMNFACT,
-     $                   MINMNUPDT
+      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP,
+     $                   KBOUND, MINMNFACT, MINMNUPDT
       REAL               HUGEVAL, TAUNAN, TEMP, TEMP2, TOL3Z
 *     ..
 *     .. External Subroutines ..
@@ -390,13 +390,13 @@
 *
       MINMNFACT = MIN( M-IOFFSET, N )
       MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
-      KMAX = MIN( KMAX, MINMNFACT )
+      KBOUND = MIN( KMAX, MINMNFACT )
       TOL3Z = SQRT( SLAMCH( 'Epsilon' ) )
       HUGEVAL = SLAMCH( 'Overflow' )
 *
 *     Compute the factorization, KK is the lomn loop index.
 *
-      DO KK = 1, KMAX
+      DO KK = 1, KBOUND
 *
          I = IOFFSET + KK
 *
@@ -674,7 +674,7 @@
 *     i.e. no condition was triggered to exit the routine.
 *     Set the number of factorized columns.
 *
-      K = KMAX
+      K = KBOUND
 *
 *     We reached the end of the loop, i.e. all KMAX columns were
 *     factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before
diff --git a/SRC/dlaqp2rk.f b/SRC/dlaqp2rk.f
index 1e3d826d7..ae2d62cac 100644
--- a/SRC/dlaqp2rk.f
+++ b/SRC/dlaqp2rk.f
@@ -355,8 +355,8 @@
       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 *     ..
 *     .. Local Scalars ..
-      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP, MINMNFACT,
-     $                   MINMNUPDT
+      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP,
+     $                   KBOUND, MINMNFACT, MINMNUPDT
       DOUBLE PRECISION   HUGEVAL, TEMP, TEMP2, TOL3Z
 *     ..
 *     .. External Subroutines ..
@@ -387,13 +387,13 @@
 *
       MINMNFACT = MIN( M-IOFFSET, N )
       MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
-      KMAX = MIN( KMAX, MINMNFACT )
+      KBOUND = MIN( KMAX, MINMNFACT )
       TOL3Z = SQRT( DLAMCH( 'Epsilon' ) )
       HUGEVAL = DLAMCH( 'Overflow' )
 *
 *     Compute the factorization, KK is the lomn loop index.
 *
-      DO KK = 1, KMAX
+      DO KK = 1, KBOUND
 *
          I = IOFFSET + KK
 *
@@ -663,7 +663,7 @@
 *     i.e. no condition was triggered to exit the routine.
 *     Set the number of factorized columns.
 *
-      K = KMAX
+      K = KBOUND
 *
 *     We reached the end of the loop, i.e. all KMAX columns were
 *     factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before
diff --git a/SRC/slaqp2rk.f b/SRC/slaqp2rk.f
index 619083358..3825e2510 100644
--- a/SRC/slaqp2rk.f
+++ b/SRC/slaqp2rk.f
@@ -355,8 +355,8 @@
       PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
 *     ..
 *     .. Local Scalars ..
-      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP, MINMNFACT,
-     $                   MINMNUPDT
+      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP,
+     $                   KBOUND, MINMNFACT, MINMNUPDT
       REAL               HUGEVAL, TEMP, TEMP2, TOL3Z
 *     ..
 *     .. External Subroutines ..
@@ -387,13 +387,13 @@
 *
       MINMNFACT = MIN( M-IOFFSET, N )
       MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
-      KMAX = MIN( KMAX, MINMNFACT )
+      KBOUND = MIN( KMAX, MINMNFACT )
       TOL3Z = SQRT( SLAMCH( 'Epsilon' ) )
       HUGEVAL = SLAMCH( 'Overflow' )
 *
 *     Compute the factorization, KK is the lomn loop index.
 *
-      DO KK = 1, KMAX
+      DO KK = 1, KBOUND
 *
          I = IOFFSET + KK
 *
@@ -663,7 +663,7 @@
 *     i.e. no condition was triggered to exit the routine.
 *     Set the number of factorized columns.
 *
-      K = KMAX
+      K = KBOUND
 *
 *     We reached the end of the loop, i.e. all KMAX columns were
 *     factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before
diff --git a/SRC/zlaqp2rk.f b/SRC/zlaqp2rk.f
index 7af2fe1f4..0e0133ecf 100644
--- a/SRC/zlaqp2rk.f
+++ b/SRC/zlaqp2rk.f
@@ -359,8 +359,8 @@
      $                   CONE = ( 1.0D+0, 0.0D+0 ) )
 *     ..
 *     .. Local Scalars ..
-      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP, MINMNFACT,
-     $                   MINMNUPDT
+      INTEGER            I, ITEMP, J, JMAXC2NRM, KK, KP,
+     $                   KBOUND, MINMNFACT, MINMNUPDT
       DOUBLE PRECISION   HUGEVAL, TAUNAN, TEMP, TEMP2, TOL3Z
 *     ..
 *     .. External Subroutines ..
@@ -391,13 +391,13 @@
 *
       MINMNFACT = MIN( M-IOFFSET, N )
       MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
-      KMAX = MIN( KMAX, MINMNFACT )
+      KBOUND = MIN( KMAX, MINMNFACT )
       TOL3Z = SQRT( DLAMCH( 'Epsilon' ) )
       HUGEVAL = DLAMCH( 'Overflow' )
 *
 *     Compute the factorization, KK is the lomn loop index.
 *
-      DO KK = 1, KMAX
+      DO KK = 1, KBOUND
 *
          I = IOFFSET + KK
 *
@@ -675,7 +675,7 @@
 *     i.e. no condition was triggered to exit the routine.
 *     Set the number of factorized columns.
 *
-      K = KMAX
+      K = KBOUND
 *
 *     We reached the end of the loop, i.e. all KMAX columns were
 *     factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before
