diff --git a/mplapack/reference/Clarnv.cpp b/mplapack/reference/Clarnv.cpp
index 5d84650a..aeb75430 100644
--- a/mplapack/reference/Clarnv.cpp
+++ b/mplapack/reference/Clarnv.cpp
@@ -36,14 +36,14 @@ void Clarnv(INTEGER const idist, INTEGER *iseed, INTEGER const n, COMPLEX *x) {
     INTEGER il = 0;
     REAL u[lv];
     INTEGER i = 0;
-    const REAL two = 2.0e+0;
+    const REAL two = 2.0;
     const REAL one = 1.0;
     const REAL zero = 0.0;
-    const REAL twopi = two * pi(zero);
+    const REAL twopi = 6.283185307179586;
     for (iv = 1; iv <= n; iv = iv + lv / 2) {
         il = min(lv / 2, n - iv + 1);
         //
-        // Call DLARUV to generate 2*IL real numbers from a uniform (0,1)
+        // Call Rlaruv to generate 2*IL real numbers from a uniform (0,1)
         // distribution (2*IL <= LV)
         //
         Rlaruv(iseed, 2 * il, u);
