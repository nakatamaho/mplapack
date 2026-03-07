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
