--- a/mplapack/reference/Rlarnv.cpp
+++ b/mplapack/reference/Rlarnv.cpp
@@ -46,7 +46,7 @@
     INTEGER i = 0;
     const REAL two = 2.0;
     const REAL one = 1.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL twopi = two * pi(one);
     for (iv = 1; iv <= n; iv = iv + lv / 2) {
         il = min(lv / 2, n - iv + 1);
         if (idist == 3) {
