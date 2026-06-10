--- Crotg.cpp_	2026-03-19 18:32:34.328480205 +0900
+++ Crotg.cpp	2026-03-19 18:34:10.119388708 +0900
@@ -34,12 +34,17 @@
 //   NAG Ltd.
 
 #include <mpblas.h>
+#include <mplapack_arithmetic_params.h>
+
+inline REAL abssq(COMPLEX ff) {
+    REAL temp;
+    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
+    return temp;
+}
 
 void Crotg(COMPLEX &a, COMPLEX const b, REAL &c, COMPLEX &s) {
     //
-    // [fable] integer, parameter :: wp = kind(1.d0)
     COMPLEX t = 0.0;
-    abssq(t) = pow2(t.real()) + pow2(t.imag());
     //
     COMPLEX f = a;
     COMPLEX g = b;
@@ -48,9 +53,10 @@
     COMPLEX r = 0.0;
     const REAL zero = 0.0;
     REAL g1 = 0.0;
-    const REAL safmax = UNHANDLED;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
     REAL rtmax = 0.0;
-    const REAL safmin = UNHANDLED;
+    const REAL safmax = ap.safmax;
     const REAL rtmin = sqrt(safmin);
     REAL g2 = 0.0;
     REAL d = 0.0;
