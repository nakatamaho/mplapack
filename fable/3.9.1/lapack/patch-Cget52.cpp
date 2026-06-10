--- Cget52.cpp~	2026-01-22 07:42:10.613767372 +0900
+++ Cget52.cpp	2026-01-22 08:11:11.421402541 +0900
@@ -91,13 +91,13 @@
     for (jvec = 1; jvec <= n; jvec = jvec + 1) {
         alphai = alpha[jvec - 1];
         betai = beta[jvec - 1];
-        abmax = max(abs1(alphai), abs1(betai));
-        if (abs1(alphai) > alfmax || abs1(betai) > betmax || abmax < one) {
+        abmax = max(cabs1(alphai), cabs1(betai));
+        if (cabs1(alphai) > alfmax || cabs1(betai) > betmax || abmax < one) {
             scale = one / max(abmax, safmin);
             alphai = scale * alphai;
             betai = scale * betai;
         }
-        scale = one / max(abs1(alphai) * bnorm, abs1(betai) * anorm, safmin);
+        scale = one / max(cabs1(alphai) * bnorm, cabs1(betai) * anorm, safmin);
         acoeff = scale * betai;
         bcoeff = scale * alphai;
         if (left) {
@@ -122,7 +122,7 @@
     for (jvec = 1; jvec <= n; jvec = jvec + 1) {
         temp1 = zero;
         for (j = 1; j <= n; j = j + 1) {
-            temp1 = max(temp1, abs1(e[(j - 1) + (jvec - 1) * lde]));
+            temp1 = max(temp1, cabs1(e[(j - 1) + (jvec - 1) * lde]));
         }
         enrmer = max(enrmer, temp1 - one);
     }
