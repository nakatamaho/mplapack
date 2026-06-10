--- a/TESTING/LIN/zchkhe_aa_2stage.f~	2021-03-26 03:25:15.000000000 +0900
+++ b/TESTING/LIN/zchkhe_aa_2stage.f	2026-01-20 12:41:38.502013178 +0900
@@ -185,7 +185,8 @@
       LOGICAL            DOTYPE( * )
       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
       COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ),
-     $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
+     $                   WORK( * ), X( * ), XACT( * )
+      DOUBLE PRECISION   RWORK( * )
 *     ..
 *
 *  =====================================================================
