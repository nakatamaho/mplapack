diff --git a/mplapack/reference/Rlarrb.cpp b/mplapack/reference/Rlarrb.cpp
index 24c24e9e..e2563386 100644
--- a/mplapack/reference/Rlarrb.cpp
+++ b/mplapack/reference/Rlarrb.cpp
@@ -66,9 +66,7 @@ void Rlarrb(INTEGER const n, REAL *d, REAL *lld, INTEGER const ifirst, INTEGER c
         return;
     }
     //
-    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + (INTEGER)2;
-    if (maxitr >= 1024)
-        maxitr = 1024; // XXX itmax can be too large for MPFR (=10^8)
+    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + 2;
     mnwdth = two * pivmin;
     //
     r = twist;
