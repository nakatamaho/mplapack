--- a/mplapack/reference/Rlassq.cpp
+++ b/mplapack/reference/Rlassq.cpp
@@ -35,6 +35,7 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
 
 void Rlassq(INTEGER const n, REAL *x, INTEGER const incx, REAL &scale, REAL &sumsq) {
     // Quick return if possible
@@ -73,10 +74,12 @@
     }
     INTEGER i = 0;
     REAL ax = 0.0;
-    const REAL tbig = 1.0 + unhandled;
-    const REAL sbig = 1.0 + unhandled;
-    const REAL tsml = 1.0 + unhandled;
-    const REAL ssml = 1.0 + unhandled;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const auto bp = mplapack::make_blue_scaling_params(ap);
+    const REAL tbig = bp.tbig;
+    const REAL sbig = bp.sbig;
+    const REAL tsml = bp.tsml;
+    const REAL ssml = bp.ssml;
     for (i = 1; i <= n; i = i + 1) {
         ax = abs(x[ix - 1]);
         if (ax > tbig) {
