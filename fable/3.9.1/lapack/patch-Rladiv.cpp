--- a/mplapack/reference/Rladiv.cpp
+++ b/mplapack/reference/Rladiv.cpp
@@ -36,11 +36,8 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-//
-//
 REAL Rladiv2(REAL const a, REAL const b, REAL const c, REAL const d, REAL const r, REAL const t) {
     REAL return_value = 0.0;
-    //
     const REAL zero = 0.0;
     REAL br = 0.0;
     if (r != zero) {
@@ -60,10 +57,7 @@
     //
 }
 
-//
-//
 void Rladiv1(REAL &a, REAL const b, REAL const c, REAL const d, REAL &p, REAL &q) {
-    //
     REAL r = d / c;
     const REAL one = 1.0;
     REAL t = one / (c + d * r);
@@ -76,7 +70,6 @@
 }
 
 void Rladiv(REAL const a, REAL const b, REAL const c, REAL const d, REAL &p, REAL &q) {
-    //
     REAL aa = a;
     REAL bb = b;
     REAL cc = c;
