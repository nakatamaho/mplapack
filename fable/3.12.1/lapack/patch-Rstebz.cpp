--- a/mplapack/reference/Rstebz.cpp
+++ b/mplapack/reference/Rstebz.cpp
@@ -230,6 +230,10 @@
         // Compute Iteration parameters
         //
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
+        if (itmax >= 100000)
+            itmax = 100000; // XXX itmax can be too large for MPFR/GMP (=10^8)
+#endif
         if (abstol <= zero) {
             atoli = ulp * tnorm;
         } else {
@@ -320,6 +324,19 @@
             //
             // Special Case -- IN=1
             //
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            if (irange == 1 || wl >= d[ibegin - 1]) {
+                nwl++;
+            }
+            if (irange == 1 || wu >= d[ibegin - 1]) {
+                nwu++;
+            }
+            if (irange == 1 || (wl < d[ibegin - 1] && wu >= d[ibegin - 1])) {
+                m++;
+                w[m - 1] = d[ibegin - 1];
+                iblock[m - 1] = jb;
+            }
+#else
             if (irange == 1 || wl >= d[ibegin - 1] - pivmin) {
                 nwl++;
             }
@@ -331,6 +348,7 @@
                 w[m - 1] = d[ibegin - 1];
                 iblock[m - 1] = jb;
             }
+#endif
         } else {
             //
             // General Case -- IN > 1
@@ -385,6 +389,10 @@
             // Compute Eigenvalues
             //
             itmax = castINTEGER((log(gu - gl + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
+            if (itmax >= 100000)
+                itmax = 100000; // XXX itmax can be too large for MPFR/GMP (=10^8)
+#endif
             Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &work[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], iout, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
             //
             // Copy Eigenvalues Into W and IBLOCK
