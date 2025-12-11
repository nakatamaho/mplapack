diff --git a/mplapack/reference/Rladiv.cpp b/mplapack/reference/Rladiv.cpp
index a8505b67..e128326e 100644
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
@@ -50,7 +55,12 @@ REAL Rladiv2(REAL const &a, REAL const &b, REAL const &c, REAL const &d, REAL co
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
@@ -63,6 +73,9 @@ void Rladiv1(REAL &a, REAL const &b, REAL const &c, REAL const &d, REAL &p, REAL
 }
 
 void Rladiv(REAL const &a, REAL const &b, REAL const &c, REAL const &d, REAL &p, REAL &q) {
+    //
+    //
+    //
     REAL aa = a;
     REAL bb = b;
     REAL cc = c;
