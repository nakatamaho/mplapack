--- a/mplapack/reference/Rlarrk.cpp
+++ b/mplapack/reference/Rlarrk.cpp
@@ -70,6 +70,10 @@
     atoli = fudge * two * pivmin;
     //
     itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined MPLAPACK_BUILD_WITH_MPFR || defined MPLAPACK_BUILD_WITH_GMP
+    if (itmax > 100000)
+        itmax = 100000; // XXX itmax can be too large for MPFR/GMP (=10^8)
+#endif
     //
     info = -1;
     //
