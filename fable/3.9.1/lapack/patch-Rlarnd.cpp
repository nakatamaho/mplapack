--- Rlarnd.cpp-	2026-01-18 08:55:00.397339978 +0900
+++ Rlarnd.cpp	2026-01-18 08:55:14.564557936 +0900
@@ -48,7 +48,7 @@
     const REAL two = 2.0;
     const REAL one = 1.0;
     REAL t2 = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL twopi = two * pi(one);
     if (idist == 1) {
         //
         // uniform (0,1)
