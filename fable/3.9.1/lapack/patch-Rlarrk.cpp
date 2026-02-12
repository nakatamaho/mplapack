--- a/mplapack/reference/Rlarrk.cpp
+++ b/mplapack/reference/Rlarrk.cpp
@@ -70,6 +70,10 @@
     atoli = fudge * two * pivmin;
     //
     itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___
+    if (itmax > 20000)
+        itmax = 20000; // XXX itmax can be too large for MPFR (=10^8)
+#endif
     //
     info = -1;
     //
