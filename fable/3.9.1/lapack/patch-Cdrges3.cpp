--- Cdrges3.cpp~	2026-01-22 15:48:55.132872152 +0900
+++ Cdrges3.cpp	2026-01-22 16:04:36.322581444 +0900
@@ -398,7 +398,7 @@
                 //
                 for (j = 1; j <= n; j = j + 1) {
                     ilabad = false;
-                    temp2 = (abs1(alpha[j - 1] - s[(j - 1) + (j - 1) * lda]) / max(safmin, abs1(alpha[j - 1]), abs1(s[(j - 1) + (j - 1) * lda])) + abs1(beta[j - 1] - t[(j - 1) + (j - 1) * lda]) / max(safmin, abs1(beta[j - 1]), abs1(t[(j - 1) + (j - 1) * lda]))) / ulp;
+                    temp2 = (cabs1(alpha[j - 1] - s[(j - 1) + (j - 1) * lda]) / max(safmin, cabs1(alpha[j - 1]), cabs1(s[(j - 1) + (j - 1) * lda])) + cabs1(beta[j - 1] - t[(j - 1) + (j - 1) * lda]) / max(safmin, cabs1(beta[j - 1]), cabs1(t[(j - 1) + (j - 1) * lda]))) / ulp;
                     //
                     if (j < n) {
                         if (s[((j + 1) - 1) + (j - 1) * lda] != zero) {
