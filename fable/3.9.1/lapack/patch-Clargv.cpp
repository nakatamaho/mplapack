diff --git a/mplapack/reference/Clargv.cpp b/mplapack/reference/Clargv.cpp
index 638a305c..02a49c5f 100644
--- a/mplapack/reference/Clargv.cpp
+++ b/mplapack/reference/Clargv.cpp
@@ -36,12 +36,6 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-inline REAL abssq(COMPLEX ff) {
-    REAL temp;
-    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
-    return temp;
-}
-
 void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, REAL *c, INTEGER const incc) {
     COMPLEX ff = 0.0;
     REAL safmin = 0.0;
@@ -73,6 +67,7 @@ void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER
     REAL dr = 0.0;
     REAL di = 0.0;
     INTEGER j = 0;
+    abs1(ff) = max(abs(ff.real()), abs(ff.imag()));
     //
     // IF( FIRST ) THEN
     // FIRST = .FALSE.
@@ -90,7 +85,7 @@ void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER
         //
         // Use identical algorithm as in Clartg
         //
-        scale = max(cabs1(f), cabs1(g));
+        scale = max(abs1(f), abs1(g));
         fs = f;
         gs = g;
         count = 0;
@@ -148,7 +143,7 @@ void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER
             cs = f2s / g2s;
             // Make sure abs(FF) = 1
             // Do complex/real division explicitly with 2 real divisions
-            if (cabs1(f) > one) {
+            if (abs1(f) > one) {
                 d = Rlapy2(f.real(), f.imag());
                 ff = COMPLEX(f.real() / d, f.imag() / d);
             } else {
