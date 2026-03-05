--- Rlatmt.cpp_	2026-01-17 18:55:18.170000865 +0900
+++ Rlatmt.cpp	2026-01-17 18:55:22.223128724 +0900
@@ -275,7 +275,8 @@
     INTEGER jku = 0;
     INTEGER jr = 0;
     REAL extra = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL two = 2.0;
+    const REAL twopi = two * pi(zero);
     REAL angle = 0.0;
     REAL c = 0.0;
     REAL s = 0.0;
