diff --git a/mplapack/reference/Cheequb.cpp b/mplapack/reference/Cheequb.cpp
index 05cb897b..16466159 100644
--- a/mplapack/reference/Cheequb.cpp
+++ b/mplapack/reference/Cheequb.cpp
@@ -41,7 +41,6 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
     bool up = false;
     const REAL zero = 0.0;
     const REAL one = 1.0;
-    const REAL two = 2.0;
     INTEGER i = 0;
     INTEGER j = 0;
     REAL tol = 0.0;
@@ -119,7 +118,7 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
         s[j - 1] = 1.0 / s[j - 1];
     }
     //
-    tol = one / sqrt(two * n);
+    tol = one / sqrt(2.0 * n);
     //
     for (iter = 1; iter <= max_iter; iter = iter + 1) {
         scale = 0.0;
@@ -149,8 +148,7 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
         // avg = s^T beta / n
         avg = 0.0;
         for (i = 1; i <= n; i = i + 1) {
-            avg += s[i - 1] * work[i - 1].real();
-            ;
+            avg += s[i - 1] * work[i - 1];
         }
         avg = avg / n;
         //
@@ -169,8 +167,8 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
             t = cabs1(a[(i - 1) + (i - 1) * lda]);
             si = s[i - 1];
             c2 = (n - 1) * t;
-            c1 = (n - 2) * (work[i - 1].real() - t * si);
-            c0 = -(t * si) * si + two * work[i - 1].real() * si - n * avg;
+            c1 = (n - 2) * (work[i - 1] - t * si);
+            c0 = -(t * si) * si + 2 * work[i - 1] * si - n * avg;
             d = c1 * c1 - 4 * c0 * c2;
             //
             if (d <= 0) {
@@ -205,7 +203,7 @@ void Cheequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
                 }
             }
             //
-            avg += (u + work[i - 1].real()) * d / n;
+            avg += (u + work[i - 1]) * d / n;
             s[i - 1] = si;
         }
     }
