--- /home/docker/Cheequb.cpp~	2026-03-21 09:44:35.275807548 +0900
+++ /home/docker/Cheequb.cpp	2026-03-21 09:45:39.972312718 +0900
@@ -41,6 +41,7 @@
     bool up = false;
     const REAL zero = 0.0;
     const REAL one = 1.0;
+    const REAL two = 2.0;
     INTEGER i = 0;
     INTEGER j = 0;
     REAL tol = 0.0;
@@ -118,7 +119,7 @@
         s[j - 1] = 1.0 / s[j - 1];
     }
     //
-    tol = one / sqrt(2.0 * n);
+    tol = one / sqrt(two * n);
     //
     for (iter = 1; iter <= max_iter; iter = iter + 1) {
         scale = 0.0;
