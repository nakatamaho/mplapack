diff --git a/mplapack/reference/Clargv.cpp b/mplapack/reference/Clargv.cpp
index 02a49c5f..e93f6bf8 100644
--- a/mplapack/reference/Clargv.cpp
+++ b/mplapack/reference/Clargv.cpp
@@ -67,7 +67,6 @@ void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER
     REAL dr = 0.0;
     REAL di = 0.0;
     INTEGER j = 0;
-    abs1(ff) = max(abs(ff.real()), abs(ff.imag()));
     //
     // IF( FIRST ) THEN
     // FIRST = .FALSE.
@@ -85,7 +84,7 @@ void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER
         //
         // Use identical algorithm as in Clartg
         //
-        scale = max(abs1(f), abs1(g));
+        scale = max(cabs1(f), cabs1(g));
         fs = f;
         gs = g;
         count = 0;
@@ -143,7 +142,7 @@ void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER
             cs = f2s / g2s;
             // Make sure abs(FF) = 1
             // Do complex/real division explicitly with 2 real divisions
-            if (abs1(f) > one) {
+            if (cabs1(f) > one) {
                 d = Rlapy2(f.real(), f.imag());
                 ff = COMPLEX(f.real() / d, f.imag() / d);
             } else {
