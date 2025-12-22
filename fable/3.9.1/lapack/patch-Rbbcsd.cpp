diff --git a/mplapack/reference/Rbbcsd.cpp b/mplapack/reference/Rbbcsd.cpp
index fe4f29be..fd4b0c3f 100644
--- a/mplapack/reference/Rbbcsd.cpp
+++ b/mplapack/reference/Rbbcsd.cpp
@@ -119,7 +119,7 @@ void Rbbcsd(const char *jobu1, const char *jobu2, const char *jobv1t, const char
     //
     INTEGER i = 0;
     const REAL zero = 0.0;
-    const REAL piover2 = 1.570796326794897;
+    const REAL piover2 = pi(zero);
     for (i = 1; i <= q; i = i + 1) {
         if (theta[i - 1] < thresh) {
             theta[i - 1] = zero;
@@ -262,7 +262,7 @@ void Rbbcsd(const char *jobu1, const char *jobu2, const char *jobv1t, const char
                 }
             } else {
                 nu = sigma21;
-                mu = sqrt(1.0 - pow2(nu));
+                mu = sqrt(one - pow2(nu));
                 if (nu < thresh) {
                     mu = one;
                     nu = zero;
