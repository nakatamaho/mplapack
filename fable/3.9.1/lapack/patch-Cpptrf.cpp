--- a/mplapack/reference/Cpptrf.cpp
+++ b/mplapack/reference/Cpptrf.cpp
@@ -82,7 +82,7 @@
             //
             // Compute U(J,J) and test for non-positive-definiteness.
             //
-            ajj = ap[jj - 1].real() - Cdotc(j - 1, &ap[jc - 1], 1, &ap[jc - 1], 1);
+            ajj = ap[jj - 1].real() - Cdotc(j - 1, &ap[jc - 1], 1, &ap[jc - 1], 1).real();
             if (ajj <= zero) {
                 ap[jj - 1] = ajj;
                 goto statement_30;
