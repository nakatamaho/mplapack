--- Clatsy.cpp
+++ Clatsy.cpp
@@ -47,8 +47,13 @@
     //
     // Initialize constants
     //
-    REAL alpha = (1.0 + sqrt(17.0)) / 8.0;
-    REAL beta = alpha - 1.0 / 1000.0;
+    REAL one = 1.0;
+    REAL two = 2.0;
+    REAL eight = 8.0;
+    REAL seventeen = 17.0;
+    REAL thousand = 1000.0;
+    REAL alpha = (one + sqrt(seventeen)) / eight;
+    REAL beta = alpha - one / thousand;
     REAL alpha3 = alpha * alpha * alpha;
     //
     // UPLO = 'U':  Upper triangular storage
@@ -61,7 +66,7 @@
     const COMPLEX eye = COMPLEX(0.0, 1.0);
     COMPLEX c = 0.0;
     COMPLEX r = 0.0;
-    if (Mlsame(uplo, "U")) {
+    if (Mlsame(uplo.elems(), "U")) {
         //
         // Fill the upper triangle of the matrix with zeros.
         //
@@ -76,7 +81,7 @@
         for (i = n; i >= n5; i = i - 5) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[(i - 1) + (i - 1) * ldx] = a;
             x[((i - 2) - 1) + (i - 1) * ldx] = b;
@@ -86,9 +91,9 @@
             x[((i - 3) - 1) + ((i - 3) - 1) * ldx] = Clarnd(2, iseed);
             x[((i - 4) - 1) + ((i - 4) - 1) * ldx] = Clarnd(2, iseed);
             if (abs(x[((i - 3) - 1) + ((i - 3) - 1) * ldx]) > abs(x[((i - 4) - 1) + ((i - 4) - 1) * ldx])) {
-                x[((i - 4) - 1) + ((i - 3) - 1) * ldx] = 2.0 * x[((i - 3) - 1) + ((i - 3) - 1) * ldx];
+                x[((i - 4) - 1) + ((i - 3) - 1) * ldx] = two * x[((i - 3) - 1) + ((i - 3) - 1) * ldx];
             } else {
-                x[((i - 4) - 1) + ((i - 3) - 1) * ldx] = 2.0 * x[((i - 4) - 1) + ((i - 4) - 1) * ldx];
+                x[((i - 4) - 1) + ((i - 3) - 1) * ldx] = two * x[((i - 4) - 1) + ((i - 4) - 1) * ldx];
             }
         }
         //
@@ -98,7 +103,7 @@
         if (i > 2) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[(i - 1) + (i - 1) * ldx] = a;
             x[((i - 2) - 1) + (i - 1) * ldx] = b;
@@ -111,9 +116,9 @@
             x[(i - 1) + (i - 1) * ldx] = Clarnd(2, iseed);
             x[((i - 1) - 1) + ((i - 1) - 1) * ldx] = Clarnd(2, iseed);
             if (abs(x[(i - 1) + (i - 1) * ldx]) > abs(x[((i - 1) - 1) + ((i - 1) - 1) * ldx])) {
-                x[((i - 1) - 1) + (i - 1) * ldx] = 2.0 * x[(i - 1) + (i - 1) * ldx];
+                x[((i - 1) - 1) + (i - 1) * ldx] = two * x[(i - 1) + (i - 1) * ldx];
             } else {
-                x[((i - 1) - 1) + (i - 1) * ldx] = 2.0 * x[((i - 1) - 1) + ((i - 1) - 1) * ldx];
+                x[((i - 1) - 1) + (i - 1) * ldx] = two * x[((i - 1) - 1) + ((i - 1) - 1) * ldx];
             }
             i = i - 2;
         } else if (i == 1) {
@@ -138,7 +143,7 @@
         for (i = 1; i <= n5; i = i + 5) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[(i - 1) + (i - 1) * ldx] = a;
             x[((i + 2) - 1) + (i - 1) * ldx] = b;
@@ -148,9 +153,9 @@
             x[((i + 3) - 1) + ((i + 3) - 1) * ldx] = Clarnd(2, iseed);
             x[((i + 4) - 1) + ((i + 4) - 1) * ldx] = Clarnd(2, iseed);
             if (abs(x[((i + 3) - 1) + ((i + 3) - 1) * ldx]) > abs(x[((i + 4) - 1) + ((i + 4) - 1) * ldx])) {
-                x[((i + 4) - 1) + ((i + 3) - 1) * ldx] = 2.0 * x[((i + 3) - 1) + ((i + 3) - 1) * ldx];
+                x[((i + 4) - 1) + ((i + 3) - 1) * ldx] = two * x[((i + 3) - 1) + ((i + 3) - 1) * ldx];
             } else {
-                x[((i + 4) - 1) + ((i + 3) - 1) * ldx] = 2.0 * x[((i + 4) - 1) + ((i + 4) - 1) * ldx];
+                x[((i + 4) - 1) + ((i + 3) - 1) * ldx] = two * x[((i + 4) - 1) + ((i + 4) - 1) * ldx];
             }
         }
         //
@@ -160,7 +165,7 @@
         if (i < n - 1) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[(i - 1) + (i - 1) * ldx] = a;
             x[((i + 2) - 1) + (i - 1) * ldx] = b;
@@ -173,9 +178,9 @@
             x[(i - 1) + (i - 1) * ldx] = Clarnd(2, iseed);
             x[((i + 1) - 1) + ((i + 1) - 1) * ldx] = Clarnd(2, iseed);
             if (abs(x[(i - 1) + (i - 1) * ldx]) > abs(x[((i + 1) - 1) + ((i + 1) - 1) * ldx])) {
-                x[((i + 1) - 1) + (i - 1) * ldx] = 2.0 * x[(i - 1) + (i - 1) * ldx];
+                x[((i + 1) - 1) + (i - 1) * ldx] = two * x[(i - 1) + (i - 1) * ldx];
             } else {
-                x[((i + 1) - 1) + (i - 1) * ldx] = 2.0 * x[((i + 1) - 1) + ((i + 1) - 1) * ldx];
+                x[((i + 1) - 1) + (i - 1) * ldx] = two * x[((i + 1) - 1) + ((i + 1) - 1) * ldx];
             }
             i += 2;
         } else if (i == n) {
