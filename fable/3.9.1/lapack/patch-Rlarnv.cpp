diff --git a/mplapack/reference/Rlarnv.cpp b/mplapack/reference/Rlarnv.cpp
index d684fdbe..c94e3ad5 100644
--- a/mplapack/reference/Rlarnv.cpp
+++ b/mplapack/reference/Rlarnv.cpp
@@ -39,7 +39,7 @@ void Rlarnv(INTEGER const idist, INTEGER *iseed, INTEGER const n, REAL *x) {
     INTEGER i = 0;
     const REAL two = 2.0;
     const REAL one = 1.0;
-    const REAL twopi = two * pi(one);
+    const REAL twopi = 6.283185307179586;
     for (iv = 1; iv <= n; iv = iv + lv / 2) {
         il = min(lv / 2, n - iv + 1);
         if (idist == 3) {
