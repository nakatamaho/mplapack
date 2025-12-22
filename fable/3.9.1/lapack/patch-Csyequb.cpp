diff --git a/mplapack/reference/Csyequb.cpp b/mplapack/reference/Csyequb.cpp
index eabaf97d..27361432 100644
--- a/mplapack/reference/Csyequb.cpp
+++ b/mplapack/reference/Csyequb.cpp
@@ -111,7 +111,7 @@ void Csyequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
         s[j - 1] = 1.0 / s[j - 1];
     }
     //
-    tol = one / sqrt(castREAL(2 * n));
+    tol = one / sqrt(2.0 * n);
     //
     for (iter = 1; iter <= max_iter; iter = iter + 1) {
         scale = 0.0;
@@ -141,7 +141,7 @@ void Csyequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
         // avg = s^T beta / n
         avg = 0.0;
         for (i = 1; i <= n; i = i + 1) {
-            avg += s[i - 1] * work[i - 1].real();
+            avg += s[i - 1] * work[i - 1];
         }
         avg = avg / n;
         //
@@ -159,9 +159,9 @@ void Csyequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
         for (i = 1; i <= n; i = i + 1) {
             t = cabs1(a[(i - 1) + (i - 1) * lda]);
             si = s[i - 1];
-            c2 = castREAL(n - 1) * t;
-            c1 = castREAL(n - 2) * (work[i - 1].real() - t * si);
-            c0 = -(t * si) * si + 2.0 * work[i - 1].real() * si - castREAL(n) * avg;
+            c2 = (n - 1) * t;
+            c1 = (n - 2) * (work[i - 1] - t * si);
+            c0 = -(t * si) * si + 2 * work[i - 1] * si - n * avg;
             d = c1 * c1 - 4 * c0 * c2;
             //
             if (d <= 0) {
@@ -196,7 +196,7 @@ void Csyequb(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, R
                 }
             }
             //
-            avg += (u + work[i - 1].real()) * d / castREA(n);
+            avg += (u + work[i - 1]) * d / n;
             s[i - 1] = si;
         }
     }
