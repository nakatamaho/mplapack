--- a/mplapack/reference/Rlarrd.cpp
+++ b/mplapack/reference/Rlarrd.cpp
@@ -201,6 +201,10 @@
         // IL through IU. The initial interval [GL,GU] from the global
         // Gerschgorin bounds GL and GU is refined by Rlaebz.
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___
+        if (itmax >= 1024)
+            itmax = 1024; // XXX itmax can be too large for MPFR (=10^8)
+#endif
         work[(n + 1) - 1] = gl;
         work[(n + 2) - 1] = gl;
         work[(n + 3) - 1] = gu;

