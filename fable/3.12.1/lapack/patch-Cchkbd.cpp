--- Cchkbd.cpp~	2026-01-22 15:48:54.830866272 +0900
+++ Cchkbd.cpp	2026-01-22 16:06:37.906154414 +0900
@@ -173,7 +173,7 @@
         m = mval[jsize - 1];
         n = nval[jsize - 1];
         mnmin = min(m, n);
-        amninv = one / max(m, n, 1);
+        amninv = one / max(m, n, (INTEGER)1);
         //
         if (nsizes != 1) {
             mtypes = min(maxtyp, ntypes);
@@ -306,7 +306,13 @@
                 //
                 // Bidiagonal, random entries
                 //
+#if defined MPLAPACK_BUILD_WITH_DD
+                temp1 = -half * log(ulp);
+#elif defined MPLAPACK_BUILD_WITH_QD
+                temp1 = -(half * half) * log(ulp);
+#else
                 temp1 = -two * log(ulp);
+#endif
                 for (j = 1; j <= mnmin; j = j + 1) {
                     bd[j - 1] = exp(temp1 * Rlarnd(2, iseed));
                     if (j < mnmin) {
@@ -514,7 +518,9 @@
             //
             for (j = 0; j <= log2ui; j = j + 1) {
+#if !defined MPLAPACK_BUILD_WITH_QD
                 Rsvdch(mnmin, bd, be, s1, temp1, iinfo);
+#endif
                 if (iinfo == 0) {
                     goto statement_140;
                 }
