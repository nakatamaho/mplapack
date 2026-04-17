--- a/mplapack/reference/Cheequb.cpp
+++ b/mplapack/reference/Cheequb.cpp
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
@@ -148,7 +149,7 @@
         // avg = s^T beta / n
         avg = 0.0;
         for (i = 1; i <= n; i = i + 1) {
-            avg += (s[i - 1] * work[i - 1]).real();
+            avg += s[i - 1] * work[i - 1].real();
         }
         avg = avg / n;
         //
@@ -168,7 +169,7 @@
             si = s[i - 1];
             c2 = (n - 1) * t;
             c1 = (n - 2) * (work[i - 1].real() - t * si);
-            c0 = -(t * si) * si + 2 * work[i - 1].real() * si - n * avg;
+            c0 = -(t * si) * si + two * work[i - 1].real() * si - n * avg;
             d = c1 * c1 - 4 * c0 * c2;
             //
             if (d <= 0) {
