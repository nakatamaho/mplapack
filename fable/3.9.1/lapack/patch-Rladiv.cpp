diff --git a/mplapack/reference/Rladiv.cpp b/mplapack/reference/Rladiv.cpp
index 69b5a6fc..c809ecee 100644
--- a/mplapack/reference/Rladiv.cpp
+++ b/mplapack/reference/Rladiv.cpp
@@ -36,8 +36,11 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+//
+//
 REAL Rladiv2(REAL const a, REAL const b, REAL const c, REAL const d, REAL const r, REAL const t) {
     REAL return_value = 0.0;
+    //
     const REAL zero = 0.0;
     REAL br = 0.0;
     if (r != zero) {
@@ -57,7 +60,10 @@ REAL Rladiv2(REAL const a, REAL const b, REAL const c, REAL const d, REAL const
     //
 }
 
+//
+//
 void Rladiv1(REAL &a, REAL const b, REAL const c, REAL const d, REAL &p, REAL &q) {
+    //
     REAL r = d / c;
     const REAL one = 1.0;
     REAL t = one / (c + d * r);
@@ -70,6 +76,7 @@ void Rladiv1(REAL &a, REAL const b, REAL const c, REAL const d, REAL &p, REAL &q
 }
 
 void Rladiv(REAL const a, REAL const b, REAL const c, REAL const d, REAL &p, REAL &q) {
+    //
     REAL aa = a;
     REAL bb = b;
     REAL cc = c;
