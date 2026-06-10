--- Clatsp.cpp
+++ Clatsp.cpp
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
     // Fill the matrix with zeros.
@@ -67,7 +72,7 @@
     const COMPLEX eye = COMPLEX(0.0, 1.0);
     COMPLEX c = 0.0;
     COMPLEX r = 0.0;
-    if (Mlsame(uplo, "U")) {
+    if (Mlsame(uplo.elems(), "U")) {
         n5 = n / 5;
         n5 = n - 5 * n5 + 1;
         //
@@ -75,7 +80,7 @@
         for (j = n; j >= n5; j = j - 5) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[jj - 1] = a;
             x[(jj - 2) - 1] = b;
@@ -89,9 +94,9 @@
             jj = jj - (j - 3);
             x[jj - 1] = Clarnd(2, iseed);
             if (abs(x[(jj + (j - 3)) - 1]) > abs(x[jj - 1])) {
-                x[(jj + (j - 4)) - 1] = 2.0 * x[(jj + (j - 3)) - 1];
+                x[(jj + (j - 4)) - 1] = two * x[(jj + (j - 3)) - 1];
             } else {
-                x[(jj + (j - 4)) - 1] = 2.0 * x[jj - 1];
+                x[(jj + (j - 4)) - 1] = two * x[jj - 1];
             }
             jj = jj - (j - 4);
         }
@@ -102,7 +107,7 @@
         if (j > 2) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[jj - 1] = a;
             x[(jj - 2) - 1] = b;
@@ -118,9 +123,9 @@
             x[jj - 1] = Clarnd(2, iseed);
             x[(jj - j) - 1] = Clarnd(2, iseed);
             if (abs(x[jj - 1]) > abs(x[(jj - j) - 1])) {
-                x[(jj - 1) - 1] = 2.0 * x[jj - 1];
+                x[(jj - 1) - 1] = two * x[jj - 1];
             } else {
-                x[(jj - 1) - 1] = 2.0 * x[(jj - j) - 1];
+                x[(jj - 1) - 1] = two * x[(jj - j) - 1];
             }
             jj = jj - j - (j - 1);
             j = j - 2;
@@ -139,7 +144,7 @@
         for (j = 1; j <= n5; j = j + 5) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[jj - 1] = a;
             x[(jj + 2) - 1] = b;
@@ -153,9 +158,9 @@
             jj += (n - j - 2);
             x[jj - 1] = Clarnd(2, iseed);
             if (abs(x[(jj - (n - j - 2)) - 1]) > abs(x[jj - 1])) {
-                x[(jj - (n - j - 2) + 1) - 1] = 2.0 * x[(jj - (n - j - 2)) - 1];
+                x[(jj - (n - j - 2) + 1) - 1] = two * x[(jj - (n - j - 2)) - 1];
             } else {
-                x[(jj - (n - j - 2) + 1) - 1] = 2.0 * x[jj - 1];
+                x[(jj - (n - j - 2) + 1) - 1] = two * x[jj - 1];
             }
             jj += (n - j - 3);
         }
@@ -166,7 +171,7 @@
         if (j < n - 1) {
             a = alpha3 * Clarnd(5, iseed);
             b = Clarnd(5, iseed) / alpha;
-            c = a - 2.0 * b * eye;
+            c = a - two * b * eye;
             r = c / beta;
             x[jj - 1] = a;
             x[(jj + 2) - 1] = b;
@@ -182,9 +187,9 @@
             x[jj - 1] = Clarnd(2, iseed);
             x[(jj + (n - j + 1)) - 1] = Clarnd(2, iseed);
             if (abs(x[jj - 1]) > abs(x[(jj + (n - j + 1)) - 1])) {
-                x[(jj + 1) - 1] = 2.0 * x[jj - 1];
+                x[(jj + 1) - 1] = two * x[jj - 1];
             } else {
-                x[(jj + 1) - 1] = 2.0 * x[(jj + (n - j + 1)) - 1];
+                x[(jj + 1) - 1] = two * x[(jj + (n - j + 1)) - 1];
             }
             jj += (n - j + 1) + (n - j);
             j += 2;
