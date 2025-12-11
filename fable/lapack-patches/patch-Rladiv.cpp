diff --git a/mplapack/reference/Rladiv.cpp b/mplapack/reference/Rladiv.cpp
index a8505b67..0f7c3b3f 100644
--- a/mplapack/reference/Rladiv.cpp
+++ b/mplapack/reference/Rladiv.cpp
@@ -29,8 +29,13 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+//
+//
 REAL Rladiv2(REAL const &a, REAL const &b, REAL const &c, REAL const &d, REAL const &r, REAL const &t) {
     REAL return_value = 0.0;
+    //
+    //
+    //
     const REAL zero = 0.0;
     REAL br = 0.0;
     if (r != zero) {
@@ -50,19 +55,27 @@ REAL Rladiv2(REAL const &a, REAL const &b, REAL const &c, REAL const &d, REAL co
     //
 }
 
+//
+//
 void Rladiv1(REAL &a, REAL const &b, REAL const &c, REAL const &d, REAL &p, REAL &q) {
+    //
+    //
+    //
     REAL r = d / c;
     const REAL one = 1.0;
     REAL t = one / (c + d * r);
-    p = Rladiv2(a, b, c, d, r, t);
+    p = dladiv2(a, b, c, d, r, t);
     a = -a;
-    q = Rladiv2(b, a, c, d, r, t);
+    q = dladiv2(b, a, c, d, r, t);
     //
     // End of Rladiv1
     //
 }
 
 void Rladiv(REAL const &a, REAL const &b, REAL const &c, REAL const &d, REAL &p, REAL &q) {
+    //
+    //
+    //
     REAL aa = a;
     REAL bb = b;
     REAL cc = c;
