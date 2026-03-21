--- a/mplapack/reference/Rlartg.cpp_	2026-03-21 10:38:17.716151903 +0900
+++ b/mplapack/reference/Rlartg.cpp	2026-03-21 10:38:53.706710284 +0900
@@ -35,12 +35,14 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
 
 void Rlartg(REAL const f, REAL const g, REAL &c, REAL &s, REAL &r) {
     const REAL one = 1.0;
-    const REAL safmin = one + unhandled;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
     REAL rtmin = sqrt(safmin);
-    const REAL safmax = one + unhandled;
+    const REAL safmax = ap.safmax;
     REAL rtmax = sqrt(safmax / 2);
     //
     REAL f1 = abs(f);
