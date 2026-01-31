--- Clattp.cpp
+++ Clattp.cpp
@@ -51,6 +51,8 @@
     REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
     REAL smlnum = unfl;
     const REAL one = 1.0;
+    const REAL half = 0.5;
+    const REAL quarter = 0.25;    
     REAL bignum = (one - ulp) / smlnum;
     Rlabad(smlnum, bignum);
     if ((imat >= 7 && imat <= 10) || imat == 18) {
@@ -228,8 +230,8 @@
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
@@ -239,7 +241,7 @@
                 work[(j + 1) - 1] = plus2;
                 work[(n + j + 1) - 1] = zero;
                 plus1 = star1 / plus2;
-                rexp = Clarnd(2, iseed);
+                rexp = Clarnd(2, iseed).real();
                 if (rexp < zero) {
                     star1 = -pow(sfac, (one - rexp)) * Clarnd(5, iseed);
                 } else {
