--- a/mplapack/reference/Rlarrd.cpp
+++ b/mplapack/reference/Rlarrd.cpp
@@ -201,6 +201,10 @@
         // IL through IU. The initial interval [GL,GU] from the global
         // Gerschgorin bounds GL and GU is refined by Rlaebz.
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined MPLAPACK_BUILD_WITH_MPFR || defined MPLAPACK_BUILD_WITH_GMP
+        if (itmax >= 100000)
+            itmax = 100000; // XXX itmax can be too large for MPFR/GMP (=10^8)
+#endif
         work[(n + 1) - 1] = gl;
         work[(n + 2) - 1] = gl;
         work[(n + 3) - 1] = gu;

