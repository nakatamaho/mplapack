--- a/mplapack/reference/Rbdsvdx.cpp
+++ b/mplapack/reference/Rbdsvdx.cpp
@@ -48,6 +48,7 @@
     //
     info = 0;
     const REAL zero = 0.0;
+    const REAL two = 2.0;
     if (!Mlsame(uplo, "U") && !lower) {
         info = -1;
     } else if (!(wantz || Mlsame(jobz, "N"))) {
@@ -107,10 +108,10 @@
         return;
     }
     //
-    REAL abstol = 2 * Rlamch("Safe Minimum");
+    REAL abstol = two * Rlamch("Safe Minimum");
     REAL ulp = Rlamch("Precision");
     REAL eps = Rlamch("Epsilon");
-    REAL sqrt2 = sqrt(2.0);
+    REAL sqrt2 = sqrt(two);
     REAL ortol = sqrt(ulp);
     //
     // Criterion for splitting is taken from Rbdsqr when singular
@@ -121,7 +122,11 @@
     const REAL ten = 10.0;
     const REAL hndrd = 100.0;
     const REAL meigth = -0.125;
+#if defined MPLAPACK_BUILD_WITH_GMP
+    REAL tol = max(ten, min(hndrd, one / sqrt(sqrt(sqrt(eps))))) * eps;
+#else
     REAL tol = max(ten, min(hndrd, pow(eps, meigth))) * eps;
+#endif
     //
     // Compute approximate maximum, minimum singular values.
     //
