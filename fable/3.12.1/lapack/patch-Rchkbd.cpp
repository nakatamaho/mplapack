--- Rchkbd.cpp~	2026-01-22 16:10:26.666902997 +0900
+++ Rchkbd.cpp	2026-01-22 16:15:01.127493207 +0900
@@ -188,7 +188,7 @@
         m = mval[jsize - 1];
         n = nval[jsize - 1];
         mnmin = min(m, n);
-        amninv = one / max(m, n, 1);
+        amninv = one / max(m, n, (INTEGER)1);
         //
         if (nsizes != 1) {
             mtypes = min(maxtyp, ntypes);
@@ -321,7 +321,11 @@
                 //
                 // Bidiagonal, random entries
                 //
+#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
+                temp1 = -half * log(ulp);
+#else
                 temp1 = -two * log(ulp);
+#endif
                 for (j = 1; j <= mnmin; j = j + 1) {
                     bd[j - 1] = exp(temp1 * Rlarnd(2, iseed));
                     if (j < mnmin) {
