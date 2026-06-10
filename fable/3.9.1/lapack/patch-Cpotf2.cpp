--- a/mplapack/reference/Cpotf2.cpp
+++ b/mplapack/reference/Cpotf2.cpp
@@ -74,7 +74,7 @@
             //
             // Compute U(J,J) and test for non-positive-definiteness.
             //
-            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1) * lda], 1, &a[(j - 1) * lda], 1);
+            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1) * lda], 1, &a[(j - 1) * lda], 1).real();
             if (ajj <= zero || Risnan(ajj)) {
                 a[(j - 1) + (j - 1) * lda] = ajj;
                 goto statement_30;
@@ -99,7 +99,7 @@
             //
             // Compute L(J,J) and test for non-positive-definiteness.
             //
-            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1)], lda, &a[(j - 1)], lda);
+            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1)], lda, &a[(j - 1)], lda).real();
             if (ajj <= zero || Risnan(ajj)) {
                 a[(j - 1) + (j - 1) * lda] = ajj;
                 goto statement_30;
