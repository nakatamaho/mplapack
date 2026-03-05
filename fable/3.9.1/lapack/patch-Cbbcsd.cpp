--- a/mplapack/reference/Cbbcsd.cpp
+++ b/mplapack/reference/Cbbcsd.cpp
@@ -126,7 +126,8 @@
     //
     INTEGER i = 0;
     const REAL zero = 0.0;
-    const REAL piover2 = 1.5707963267948966192313216916397514421;
+    const REAL two = 2.0;
+    const REAL piover2 = pi(zero) / two;
     for (i = 1; i <= q; i = i + 1) {
         if (theta[i - 1] < thresh) {
             theta[i - 1] = zero;
@@ -269,7 +270,7 @@
                 }
             } else {
                 nu = sigma21;
-                mu = sqrt(1.0 - pow2(nu));
+                mu = sqrt(one - pow2(nu));
                 if (nu < thresh) {
                     mu = one;
                     nu = zero;
