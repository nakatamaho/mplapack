--- a/mplapack/reference/Rlarrb.cpp
+++ b/mplapack/reference/Rlarrb.cpp
@@ -73,7 +73,11 @@
         return;
     }
     //
-    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + 2;
+    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + (INTEGER)2;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
+    if (maxitr >= 100000)
+        maxitr = 100000; // XXX itmax can be too large for MPFR/GMP (=10^8)
+#endif
     mnwdth = two * pivmin;
     //
     r = twist;
