diff --git a/mplapack/reference/Cbbcsd.cpp b/mplapack/reference/Cbbcsd.cpp
index 9e133841..b152bb91 100644
--- a/mplapack/reference/Cbbcsd.cpp
+++ b/mplapack/reference/Cbbcsd.cpp
@@ -120,7 +120,7 @@ void Cbbcsd(const char *jobu1, const char *jobu2, const char *jobv1t, const char
     //
     INTEGER i = 0;
     const REAL zero = 0.0;
-    const REAL piover2 = pi(zero) / 2.0;
+    const REAL piover2 = 1.570796326794897;
     for (i = 1; i <= q; i = i + 1) {
         if (theta[i - 1] < thresh) {
             theta[i - 1] = zero;
@@ -263,7 +263,7 @@ void Cbbcsd(const char *jobu1, const char *jobu2, const char *jobv1t, const char
                 }
             } else {
                 nu = sigma21;
-                mu = sqrt(one - pow2(nu));
+                mu = sqrt(1.0 - pow2(nu));
                 if (nu < thresh) {
                     mu = one;
                     nu = zero;
