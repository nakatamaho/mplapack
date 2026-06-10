--- Rnrm2.cpp_	2026-03-19 18:40:54.918591972 +0900
+++ Rnrm2.cpp	2026-03-19 18:45:16.912866364 +0900
@@ -34,6 +34,7 @@
 //   NAG Ltd.
 
 #include <mpblas.h>
+#include <mplapack_arithmetic_params.h>
 
 REAL Rnrm2(INTEGER const n, REAL *x, INTEGER const incx) {
     REAL return_value = 0.0;
@@ -67,10 +68,12 @@
     }
     INTEGER i = 0;
     REAL ax = 0.0;
-    const REAL tbig = UNHANDLED;
-    const REAL sbig = UNHANDLED;
-    const REAL tsml = UNHANDLED;
-    const REAL ssml = UNHANDLED;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const auto bp = mplapack::make_blue_scaling_params(ap);
+    const REAL tbig = bp.tbig;
+    const REAL sbig = bp.sbig;
+    const REAL tsml = bp.tsml;
+    const REAL ssml = bp.ssml;
     for (i = 1; i <= n; i = i + 1) {
         ax = abs(x[ix - 1]);
         if (ax > tbig) {
@@ -89,7 +92,7 @@
     // Combine abig and amed or amed and asml if more than one
     // accumulator was used.
     //
-    const REAL maxn = Rlamch("O");
+    const REAL maxn = ap.rmax;
     REAL ymin = 0.0;
     REAL ymax = 0.0;
     if (abig > zero) {
