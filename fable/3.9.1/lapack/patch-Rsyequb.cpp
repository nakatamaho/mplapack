--- a/mplapack/reference/Rsyequb.cpp
+++ b/mplapack/reference/Rsyequb.cpp
@@ -40,6 +40,7 @@
     bool up = false;
     const REAL zero = 0.0;
     const REAL one = 1.0;
+    const REAL two = 2.0;
     INTEGER i = 0;
     INTEGER j = 0;
     REAL tol = 0.0;
@@ -117,7 +118,7 @@
         s[j - 1] = 1.0 / s[j - 1];
     }
     //
-    tol = one / sqrt(2.0 * n);
+    tol = one / sqrt(two * n);
     //
     for (iter = 1; iter <= max_iter; iter = iter + 1) {
         scale = 0.0;
