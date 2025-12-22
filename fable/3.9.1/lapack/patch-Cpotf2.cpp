diff --git a/mplapack/reference/Cpotf2.cpp b/mplapack/reference/Cpotf2.cpp
index 5243f821..13500383 100644
--- a/mplapack/reference/Cpotf2.cpp
+++ b/mplapack/reference/Cpotf2.cpp
@@ -67,7 +67,7 @@ void Cpotf2(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, IN
             //
             // Compute U(J,J) and test for non-positive-definiteness.
             //
-            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1) * lda], 1, &a[(j - 1) * lda], 1).real();
+            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1) * lda], 1, &a[(j - 1) * lda], 1);
             if (ajj <= zero || Risnan(ajj)) {
                 a[(j - 1) + (j - 1) * lda] = ajj;
                 goto statement_30;
@@ -92,7 +92,7 @@ void Cpotf2(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, IN
             //
             // Compute L(J,J) and test for non-positive-definiteness.
             //
-            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1)], lda, &a[(j - 1)], lda).real();
+            ajj = a[(j - 1) + (j - 1) * lda].real() - Cdotc(j - 1, &a[(j - 1)], lda, &a[(j - 1)], lda);
             if (ajj <= zero || Risnan(ajj)) {
                 a[(j - 1) + (j - 1) * lda] = ajj;
                 goto statement_30;
