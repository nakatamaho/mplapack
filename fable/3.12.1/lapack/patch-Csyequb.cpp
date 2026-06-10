--- a/mplapack/reference/Csyequb.cpp
+++ b/mplapack/reference/Csyequb.cpp
@@ -41,6 +41,7 @@
     bool up = false;
     const REAL zero = 0.0;
     const REAL one = 1.0;
+    const REAL two = 2.0;
     INTEGER i = 0;
     INTEGER j = 0;
     REAL tol = 0.0;
@@ -118,7 +119,7 @@
         s[j - 1] = 1.0 / s[j - 1];
     }
     //
-    tol = one / sqrt(2.0 * n);
+    tol = one / sqrt(two * n);
     //
     for (iter = 1; iter <= max_iter; iter = iter + 1) {
         scale = 0.0;
@@ -166,9 +167,9 @@
         for (i = 1; i <= n; i = i + 1) {
             t = cabs1(a[(i - 1) + (i - 1) * lda]);
             si = s[i - 1];
-            c2 = (n - 1) * t;
-            c1 = (n - 2) * (work[i - 1].real() - t * si);
-            c0 = -(t * si) * si + 2 * work[i - 1].real() * si - n * avg;
+            c2 = castREAL(n - 1) * t;
+            c1 = castREAL(n - 2) * (work[i - 1].real() - t * si);
+            c0 = -(t * si) * si + two * work[i - 1].real() * si - castREAL(n) * avg;
             d = c1 * c1 - 4 * c0 * c2;
             //
             if (d <= 0) {
@@ -203,7 +204,7 @@
                 }
             }
             //
-            avg += (u + work[i - 1].real()) * d / n;
+            avg += (u + work[i - 1].real()) * d / castREAL(n);
             s[i - 1] = si;
         }
     }
