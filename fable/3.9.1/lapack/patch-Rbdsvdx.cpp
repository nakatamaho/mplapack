diff --git a/mplapack/reference/Rbdsvdx.cpp b/mplapack/reference/Rbdsvdx.cpp
index c7b828ba..a15f93a5 100644
--- a/mplapack/reference/Rbdsvdx.cpp
+++ b/mplapack/reference/Rbdsvdx.cpp
@@ -48,6 +48,7 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
     //
     info = 0;
     const REAL zero = 0.0;
+    const REAL two = 2.0;
     if (!Mlsame(uplo, "U") && !lower) {
         info = -1;
     } else if (!(wantz || Mlsame(jobz, "N"))) {
@@ -107,10 +108,10 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
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
