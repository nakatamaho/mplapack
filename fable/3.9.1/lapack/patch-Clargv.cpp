--- a/mplapack/reference/Clargv.cpp
+++ b/mplapack/reference/Clargv.cpp
@@ -36,6 +36,12 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL abssq(COMPLEX ff) {
+    REAL temp;
+    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
+    return temp;
+}
+
 void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, REAL *c, INTEGER const incc) {
     COMPLEX ff = 0.0;
     REAL safmin = 0.0;
@@ -67,7 +73,6 @@
     REAL dr = 0.0;
     REAL di = 0.0;
     INTEGER j = 0;
-    abs1(ff) = max(abs(ff.real()), abs(ff.imag()));
     //
     // IF( FIRST ) THEN
     // FIRST = .FALSE.
@@ -85,7 +90,7 @@
         //
         // Use identical algorithm as in Clartg
         //
-        scale = max(abs1(f), abs1(g));
+        scale = max(cabs1(f), cabs1(g));
         fs = f;
         gs = g;
         count = 0;
@@ -143,7 +148,7 @@
             cs = f2s / g2s;
             // Make sure abs(FF) = 1
             // Do complex/real division explicitly with 2 real divisions
-            if (abs1(f) > one) {
+            if (cabs1(f) > one) {
                 d = Rlapy2(f.real(), f.imag());
                 ff = COMPLEX(f.real() / d, f.imag() / d);
             } else {
