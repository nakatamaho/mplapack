--- a/mplapack/reference/Rbbcsd.cpp
+++ b/mplapack/reference/Rbbcsd.cpp
@@ -126,7 +130,7 @@
     //
     INTEGER i = 0;
     const REAL zero = 0.0;
-    const REAL piover2 = 1.5707963267948966192313216916397514421;
+    const REAL piover2 = pi(zero) / 2.0;
     for (i = 1; i <= q; i = i + 1) {
         if (theta[i - 1] < thresh) {
             theta[i - 1] = zero;
@@ -269,7 +273,7 @@
                 }
             } else {
                 nu = sigma21;
-                mu = sqrt(1.0 - pow2(nu));
+                mu = sqrt(one - pow2(nu));
                 if (nu < thresh) {
                     mu = one;
                     nu = zero;
