--- a/mplapack/reference/Clartg.cpp	2026-03-21 11:38:06.661888232 +0900
+++ b/mplapack/reference/Clartg.cpp	2026-03-21 11:42:48.544846415 +0900
@@ -35,17 +35,19 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
 
 void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
     COMPLEX t = 0.0;
-    const REAL safmin = 1.0 + unhandled;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
     REAL rtmin = sqrt(safmin);
     //
     const COMPLEX czero = COMPLEX(0.0, 0.0);
     const REAL one = 1.0;
     const REAL zero = 0.0;
     REAL g1 = 0.0;
-    const REAL safmax = 1.0 + unhandled;
+    const REAL safmax = ap.safmax;
     REAL rtmax = 0.0;
     REAL g2 = 0.0;
     REAL d = 0.0;
