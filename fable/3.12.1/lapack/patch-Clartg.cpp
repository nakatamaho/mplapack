--- a/mplapack/reference/Clartg.cpp_	2026-03-21 10:37:02.941040414 +0900
+++ b/mplapack/reference/Clartg.cpp	2026-03-21 10:37:12.767182475 +0900
@@ -35,17 +35,19 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
 
 void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
     COMPLEX t = 0.0;
     const REAL one = 1.0;
-    const REAL safmin = one + unhandled;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
     REAL rtmin = sqrt(safmin);
     //
     REAL czero = 0.0;
     const REAL zero = 0.0;
     REAL g1 = 0.0;
-    const REAL safmax = one + unhandled;
+    const REAL safmax = ap.safmax;
     REAL rtmax = 0.0;
     REAL g2 = 0.0;
     REAL d = 0.0;
