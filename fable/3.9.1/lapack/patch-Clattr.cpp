--- Clattr.cpp
+++ Clattr.cpp
@@ -51,6 +51,8 @@
     REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
     REAL smlnum = unfl;
     const REAL one = 1.0;
+    const REAL half = 0.5;
+    const REAL quarter = 0.25;
     REAL bignum = (one - ulp) / smlnum;
     Rlabad(smlnum, bignum);
     if ((imat >= 7 && imat <= 10) || imat == 18) {
@@ -213,8 +215,8 @@
         //
         // where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).
         //
-        star1 = 0.25 * Clarnd(5, iseed);
-        sfac = 0.5;
+        star1 = quarter * Clarnd(5, iseed);
+        sfac = half;
         plus1 = sfac * Clarnd(5, iseed);
         for (j = 1; j <= n; j = j + 2) {
             plus2 = star1 / plus1;
@@ -235,7 +237,7 @@
         //
         x = sqrt(cndnum) - 1 / sqrt(cndnum);
         if (n > 2) {
-            y = sqrt(2.0 / (n - 2)) * x;
+            y = sqrt(two / REAL(static_cast<double>(n - 2))) * x; //DD, QD don't have constructor dd_real(long int) and doesn't harm for other environment
         } else {
             y = zero;
         }
@@ -272,7 +274,7 @@
         if (upper) {
             for (j = 1; j <= n - 1; j = j + 1) {
                 ra = a[(j - 1) + ((j + 1) - 1) * lda];
-                rb = 2.0;
+                rb = two;
                 Crotg(ra, rb, c, s);
                 //
                 // Multiply by [ c  s; -conjg(s)  c] on the left.
@@ -294,7 +296,7 @@
         } else {
             for (j = 1; j <= n - 1; j = j + 1) {
                 ra = a[((j + 1) - 1) + (j - 1) * lda];
-                rb = 2.0;
+                rb = two;
                 Crotg(ra, rb, c, s);
                 s = conj(s);
                 //
@@ -538,7 +540,7 @@
                 a[((j - 1) - 1) * lda] = -(tscal / castREAL(n + 1)) / castREAL(n + 2);
                 a[((j - 1) - 1) + ((j - 1) - 1) * lda] = one;
                 b[(j - 1) - 1] = texp * castREAL(n * n + n - 1);
-                texp = texp * 2.0;
+                texp = texp * two;
             }
             b[1 - 1] = (castREAL(n + 1) / castREAL(n + 2)) * tscal;
         } else {
@@ -549,7 +551,7 @@
                 a[(n - 1) + ((j + 1) - 1) * lda] = -(tscal / castREAL(n + 1)) / castREAL(n + 2);
                 a[((j + 1) - 1) + ((j + 1) - 1) * lda] = one;
                 b[(j + 1) - 1] = texp * castREAL(n * n + n - 1);
-                texp = texp * 2.0;
+                texp = texp * two;
             }
             b[n - 1] = (castREAL(n + 1) / castREAL(n + 2)) * tscal;
         }
