--- a/mplapack/reference/Rstebz.cpp	2026-01-21 21:03:48.611266540 +0900
+++ b/mplapack/reference/Rstebz.cpp	2026-01-21 21:03:48.618266651 +0900
@@ -55,6 +55,7 @@
     REAL gl = 0.0;
     REAL tmp2 = 0.0;
     REAL tnorm = 0.0;
+    REAL default_atol;
     const REAL fudge = 2.1;
     const REAL two = 2.0;
     INTEGER itmax = 0;
@@ -230,10 +231,11 @@
         // Compute Iteration parameters
         //
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+        default_atol = ulp * tnorm;
         if (abstol <= zero) {
             atoli = ulp * tnorm;
         } else {
-            atoli = abstol;
+            atoli = max(abstol, default_atol);
         }
         //
         work[(n + 1) - 1] = gl;
@@ -281,10 +283,11 @@
             tnorm = max(tnorm, abs(d[j - 1]) + abs(e[(j - 1) - 1]) + abs(e[j - 1]));
         }
         //
+        default_atol = ulp * tnorm;
         if (abstol <= zero) {
             atoli = ulp * tnorm;
         } else {
-            atoli = abstol;
+            atoli = max(abstol, default_atol);
         }
         //
         if (irange == 2) {
@@ -353,10 +356,11 @@
             //
             // Compute ATOLI for the current submatrix
             //
+            default_atol = ulp * max(abs(gl), abs(gu));
             if (abstol <= zero) {
-                atoli = ulp * max(abs(gl), abs(gu));
+                atoli = default_atol;
             } else {
-                atoli = abstol;
+                atoli = max(abstol, default_atol);
             }
             //
             if (irange > 1) {
