diff --git a/mplapack/reference/Clargv.cpp b/mplapack/reference/Clargv.cpp
index e93f6bf8..638a305c 100644
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
