--- Clatmt.cpp
+++ Clatmt.cpp
@@ -283,7 +283,8 @@ void Clatmt(INTEGER const m, INTEGER const n, const char *dist, INTEGER *iseed,
     INTEGER jku = 0;
     INTEGER jr = 0;
     COMPLEX extra = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL two = 2.0;
+    const REAL twopi = two * pi(zero);
     REAL angle = 0.0;
     COMPLEX c = 0.0;
     COMPLEX s = 0.0;
