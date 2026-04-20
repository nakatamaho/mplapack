--- Clattp.cpp	2026-03-23 17:51:50.188955148 +0900
+++ Clattp.cpp_	2026-03-23 17:50:58.181143181 +0900
@@ -51,6 +51,8 @@
     REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
     REAL smlnum = unfl;
     const REAL one = 1.0;
+    const REAL half = 0.5;
+    const REAL quarter = 0.25;
     REAL bignum = (one - ulp) / smlnum;
     if ((imat >= 7 && imat <= 10) || imat == 18) {
         diag = "U";
@@ -229,8 +227,8 @@
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
