diff --git a/mplapack/reference/Rlarrk.cpp b/mplapack/reference/Rlarrk.cpp
index 07f04fd0..fe6fc315 100644
--- a/mplapack/reference/Rlarrk.cpp
+++ b/mplapack/reference/Rlarrk.cpp
@@ -63,8 +63,6 @@ void Rlarrk(INTEGER const n, INTEGER const iw, REAL const &gl, REAL const &gu, R
     atoli = fudge * two * pivmin;
     //
     itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
-    if (itmax > 1024)
-        itmax = 1024; // XXX itmax can be too large for MPFR (=10^8)
     //
     info = -1;
     //
