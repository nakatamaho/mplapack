--- a/mplapack/reference/Rlartg.cpp
+++ b/mplapack/reference/Rlartg.cpp
@@ -35,11 +35,13 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
 
 void Rlartg(REAL const f, REAL const g, REAL &c, REAL &s, REAL &r) {
-    const REAL safmin = 1.0 + unhandled;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
     REAL rtmin = sqrt(safmin);
-    const REAL safmax = 1.0 + unhandled;
+    const REAL safmax = ap.safmax;
     REAL rtmax = sqrt(safmax / 2);
     //
     REAL f1 = abs(f);
