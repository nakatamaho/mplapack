--- a/mplapack/reference/Clargv.cpp	2026-03-28 17:15:39.952228102 +0900
+++ b/mplapack/reference/Clargv.cpp	2026-03-28 17:01:08.929227721 +0900
@@ -36,6 +36,18 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL abssq(COMPLEX ff) {
+    REAL temp;
+    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
+    return temp;
+}
+
+inline REAL cabsmax(COMPLEX z) {
+    REAL temp;
+    temp = max(abs(z.real()), abs(z.imag()));
+    return temp;
+}
+
 void Clargv(INTEGER const n, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, REAL *c, INTEGER const incc) {
     COMPLEX ff = 0.0;
     REAL safmin = 0.0;
@@ -67,7 +79,6 @@
     REAL dr = 0.0;
     REAL di = 0.0;
     INTEGER j = 0;
-    cabsmax(ff) = max(abs(ff.real()), abs(ff.imag()));
     //
     // IF( FIRST ) THEN
     // FIRST = .FALSE.
