--- a/mplapack/reference/Rlacon.cpp_	2026-01-21 21:03:48.611266540 +0900
+++ b/mplapack/reference/Rlacon.cpp	2026-01-21 21:03:48.618266651 +0900
+++ Rlacon.cpp	2026-01-21 21:03:48.634266905 +0900
@@ -37,14 +37,14 @@
 #include <mplapack.h>
 
 void Rlacon(INTEGER const n, REAL *v, REAL *x, INTEGER *isgn, REAL &est, INTEGER &kase) {
-    REAL altsgn = 0.0;
-    REAL estold = 0.0;
-    INTEGER i = 0;
-    INTEGER iter = 0;
-    INTEGER j = 0;
-    INTEGER jlast = 0;
-    INTEGER jump = 0;
-    REAL temp = 0.0;
+    static REAL altsgn = 0.0;
+    static REAL estold = 0.0;
+    static INTEGER i = 0;
+    static INTEGER iter = 0;
+    static INTEGER j = 0;
+    static INTEGER jlast = 0;
+    static INTEGER jump = 0;
+    static REAL temp = 0.0;
     const REAL one = 1.0;
     const REAL zero = 0.0;
     const INTEGER itmax = 5;
