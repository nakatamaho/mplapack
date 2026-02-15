--- a/mplapack/reference/Rstebz.cpp
+++ b/mplapack/reference/Rstebz.cpp
@@ -1,5 +1,5 @@
 /*
- * Copyright (c) 2008-2025
+ * Copyright (c) 2008-2026
  *      Nakata, Maho
  *      All rights reserved.
  *
@@ -86,6 +86,7 @@
     INTEGER iw = 0;
     INTEGER ie = 0;
     INTEGER itmp1 = 0;
+    REAL default_atol;
     //
     info = 0;
     //
@@ -230,11 +231,32 @@
         // Compute Iteration parameters
         //
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___
+        if (itmax >= 20000)
+            itmax = 20000; // XXX itmax can be too large for MPFR (=10^8)
+
+        // Choose ATOLI. In wide-exponent MPFR profiles, extremely tiny positive ABSTOL
+        // (e.g. ~ UNFL) can cause excessive bisection iterations. Clamp it to a sane floor.
+        const bool mpfr_wide_exp = (mpfr_get_emax() - mpfr_get_emin() >= 20000);
+
+        default_atol = ulp * tnorm;
+        if (mpfr_wide_exp && abstol > zero) {
+            atoli = max(abstol, default_atol);
+        } else {
+            if (abstol <= zero) {
+                atoli = default_atol;
+            } else {
+                atoli = abstol;
+            }
+        }
+#else
         if (abstol <= zero) {
             atoli = ulp * tnorm;
         } else {
             atoli = abstol;
         }
+#endif
+
         //
         work[(n + 1) - 1] = gl;
         work[(n + 2) - 1] = gl;
@@ -281,11 +303,28 @@
             tnorm = max(tnorm, abs(d[j - 1]) + abs(e[(j - 1) - 1]) + abs(e[j - 1]));
         }
         //
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___
+        // Choose ATOLI. In wide-exponent MPFR profiles, extremely tiny positive ABSTOL
+        // (e.g. ~ UNFL) can cause excessive bisection iterations. Clamp it to a sane floor.
+        const bool mpfr_wide_exp = (mpfr_get_emax() - mpfr_get_emin() >= 20000);
+
+        default_atol = ulp * tnorm;
+        if (mpfr_wide_exp && abstol > zero) {
+            atoli = max(abstol, default_atol);
+        } else {
+            if (abstol <= zero) {
+                atoli = default_atol;
+            } else {
+                atoli = abstol;
+            }
+        }
+#else
         if (abstol <= zero) {
             atoli = ulp * tnorm;
         } else {
             atoli = abstol;
         }
+#endif
         //
         if (irange == 2) {
             wl = vl;
@@ -353,11 +392,28 @@
             //
             // Compute ATOLI for the current submatrix
             //
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___
+            // Choose ATOLI. In wide-exponent MPFR profiles, extremely tiny positive ABSTOL
+            // (e.g. ~ UNFL) can cause excessive bisection iterations. Clamp it to a sane floor.
+            const bool mpfr_wide_exp = (mpfr_get_emax() - mpfr_get_emin() >= 20000);
+
+            default_atol = ulp * max(abs(gl), abs(gu));
+            if (mpfr_wide_exp && abstol > zero) {
+                atoli = max(abstol, default_atol);
+            } else {
+                if (abstol <= zero) {
+                    atoli = default_atol;
+                } else {
+                    atoli = abstol;
+                }
+            }
+#else
             if (abstol <= zero) {
                 atoli = ulp * max(abs(gl), abs(gu));
             } else {
                 atoli = abstol;
             }
+#endif
             //
             if (irange > 1) {
                 if (gu < wl) {
@@ -385,6 +441,10 @@
             // Compute Eigenvalues
             //
             itmax = castINTEGER((log(gu - gl + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___
+            if (itmax >= 20000)
+                itmax = 20000; // XXX itmax can be too large for MPFR (=10^8)
+#endif
             Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &work[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], iout, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
             //
             // Copy Eigenvalues Into W and IBLOCK
