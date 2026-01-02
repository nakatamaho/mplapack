--- a/mplapack/reference/Rlarrk.cpp
+++ b/mplapack/reference/Rlarrk.cpp
@@ -70,6 +70,8 @@
     atoli = fudge * two * pivmin;
     //
     itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
+    if (itmax > 1024)
+        itmax = 1024; // XXX itmax can be too large for MPFR (=10^8)
     //
     info = -1;
     //
