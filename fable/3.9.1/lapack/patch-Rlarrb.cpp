--- a/mplapack/reference/Rlarrb.cpp
+++ b/mplapack/reference/Rlarrb.cpp
@@ -73,7 +73,9 @@
         return;
     }
     //
-    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + 2;
+    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + (INTEGER)2;
+    if (maxitr >= 1024)
+        maxitr = 1024; // XXX itmax can be too large for MPFR (=10^8)
     mnwdth = two * pivmin;
     //
     r = twist;
