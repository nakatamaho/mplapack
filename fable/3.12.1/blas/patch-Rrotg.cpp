--- Rrotg.cpp_	2026-03-19 18:35:12.807645987 +0900
+++ Rrotg.cpp	2026-03-19 18:45:50.288921215 +0900
@@ -34,16 +34,17 @@
 //   NAG Ltd.
 
 #include <mpblas.h>
+#include <mplapack_arithmetic_params.h>
 
 void Rrotg(REAL &a, REAL &b, REAL &c, REAL &s) {
     //
-    // [fable] integer, parameter :: wp = kind(1.d0)
     REAL anorm = abs(a);
     REAL bnorm = abs(b);
     const REAL zero = 0.0;
     const REAL one = 1.0;
-    const REAL safmax = UNHANDLED;
-    const REAL safmin = UNHANDLED;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
+    const REAL safmax = ap.safmax;
     REAL scl = 0.0;
     REAL sigma = 0.0;
     REAL r = 0.0;
