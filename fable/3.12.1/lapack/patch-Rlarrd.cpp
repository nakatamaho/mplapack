--- a/mplapack/reference/Rlarrd.cpp
+++ b/mplapack/reference/Rlarrd.cpp
@@ -195,6 +195,10 @@
         // IL through IU. The initial interval [GL,GU] from the global
         // Gerschgorin bounds GL and GU is refined by Rlaebz.
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
+        if (itmax >= 100000)
+            itmax = 100000; // XXX itmax can be too large for MPFR/GMP (=10^8)
+#endif
         work[(n + 1) - 1] = gl;
         work[(n + 2) - 1] = gl;
         work[(n + 3) - 1] = gu;
