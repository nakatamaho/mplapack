diff --git a/mplapack/reference/Cheequb.cpp b/mplapack/reference/Cheequb.cpp
index 16466159..0a17f744 100644
--- a/mplapack/reference/Cheequb.cpp
+++ b/mplapack/reference/Cheequb.cpp
@@ -41,6 +41,7 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
     bool up = false;
     const REAL zero = 0.0;
     const REAL one = 1.0;
+    const REAL two = 2.0;
     INTEGER i = 0;
     INTEGER j = 0;
     REAL tol = 0.0;
@@ -118,7 +119,7 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
         s[j - 1] = 1.0 / s[j - 1];
     }
     //
-    tol = one / sqrt(2.0 * n);
+    tol = one / sqrt(two * castREAL(n));
     //
     for (iter = 1; iter <= max_iter; iter = iter + 1) {
         scale = 0.0;
