--- Clatms.cpp-	2026-01-18 08:46:45.100020982 +0900
+++ Clatms.cpp	2026-01-18 08:47:20.670520119 +0900
@@ -145,7 +145,7 @@
             givens = true;
         }
     } else {
-        if (2 * llb < m) {
+        if ((INTEGER)2 * llb < m) {
             givens = true;
         }
     }
@@ -283,7 +279,8 @@
     INTEGER jku = 0;
     INTEGER jr = 0;
     COMPLEX extra = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL two = 2.0;
+    const REAL twopi = two * pi(temp);
     REAL angle = 0.0;
     COMPLEX c = 0.0;
     COMPLEX s = 0.0;
