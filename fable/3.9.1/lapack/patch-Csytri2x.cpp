--- a/mplapack/reference/Csytri2x.cpp
+++ b/mplapack/reference/Csytri2x.cpp
@@ -130,7 +130,7 @@
         while (k <= n) {
             if (ipiv[k - 1] > 0) {
                 // 1 x 1 diagonal NNB
-                work[(k - 1) + (invd - 1) * ldwork] = 1 / a[(k - 1) + (k - 1) * lda];
+                work[(k - 1) + (invd - 1) * ldwork] = one / a[(k - 1) + (k - 1) * lda];
                 work[(k - 1) + ((invd + 1) - 1) * ldwork] = 0.0;
                 k++;
             } else {
@@ -309,7 +309,7 @@
         while (k >= 1) {
             if (ipiv[k - 1] > 0) {
                 // 1 x 1 diagonal NNB
-                work[(k - 1) + (invd - 1) * ldwork] = 1 / a[(k - 1) + (k - 1) * lda];
+                work[(k - 1) + (invd - 1) * ldwork] = one / a[(k - 1) + (k - 1) * lda];
                 work[(k - 1) + ((invd + 1) - 1) * ldwork] = 0.0;
                 k = k - 1;
             } else {
