diff --git a/mplapack/reference/Cpptrf.cpp b/mplapack/reference/Cpptrf.cpp
index 44e9d625..3543d4a2 100644
--- a/mplapack/reference/Cpptrf.cpp
+++ b/mplapack/reference/Cpptrf.cpp
@@ -75,7 +75,7 @@ void Cpptrf(const char *uplo, INTEGER const n, COMPLEX *ap, INTEGER &info) {
             //
             // Compute U(J,J) and test for non-positive-definiteness.
             //
-            ajj = ap[jj - 1].real() - Cdotc(j - 1, &ap[jc - 1], 1, &ap[jc - 1], 1).real();
+            ajj = ap[jj - 1].real() - Cdotc(j - 1, &ap[jc - 1], 1, &ap[jc - 1], 1);
             if (ajj <= zero) {
                 ap[jj - 1] = ajj;
                 goto statement_30;
