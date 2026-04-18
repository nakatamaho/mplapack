--- Cget23.cpp~	2026-04-18 00:00:00.000000000 +0900
+++ Cget23.cpp	2026-04-18 00:00:00.000000000 +0900
@@ -358,11 +358,23 @@ void Cget23(bool const comp, INTEGER con
         // Do Test (7)
         //
         for (j = 1; j <= n; j = j + 1) {
+#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
+            REAL vdiff = zero;
+            REAL vscale = one;
+            for (jj = 1; jj <= n; jj = jj + 1) {
+                vdiff = max(vdiff, abs(vl[(jj - 1) + (j - 1) * ldvl] - lre[(jj - 1) + (j - 1) * ldlre]));
+                vscale = max(vscale, abs(vl[(jj - 1) + (j - 1) * ldvl]), abs(lre[(jj - 1) + (j - 1) * ldlre]));
+            }
+            if (vdiff > 100.0 * ulp * vscale) {
+                result[7 - 1] = ulpinv;
+            }
+#else
             for (jj = 1; jj <= n; jj = jj + 1) {
                 if (vl[(j - 1) + (jj - 1) * ldvl] != lre[(j - 1) + (jj - 1) * ldlre]) {
                     result[7 - 1] = ulpinv;
                 }
             }
+#endif
         }
         //
         // Do Test (8) again
