diff --git a/mplapack/test/eig/common/Rcsdts.cpp b/mplapack/test/eig/common/Rcsdts.cpp
index 13160d44..2caca622 100644
--- Rcsdts.cpp
+++ Rcsdts.cpp
@@ -184,8 +184,9 @@ void Rcsdts(INTEGER const m, INTEGER const p, INTEGER const q, REAL *x, REAL *xf
     // Check sorting
     //
     const REAL realzero = 0.0;
+    const REAL two = 2.0;
     result[9 - 1] = realzero;
-    const REAL piover2 = 1.5707963267948966192313216916397514421;
+    const REAL piover2 = pi(realzero) / two;
     for (i = 1; i <= r; i = i + 1) {
         if (theta[i - 1] < realzero || theta[i - 1] > piover2) {
             result[9 - 1] = ulpinv;
