--- a/mplapack/reference/Clacon.cpp_	2026-01-21 21:03:48.611266540 +0900
+++ b/mplapack/reference/Clacon.cpp	2026-01-21 21:03:48.618266651 +0900
@@ -37,16 +37,16 @@
 #include <mplapack.h>
 
 void Clacon(INTEGER const n, COMPLEX *v, COMPLEX *x, REAL &est, INTEGER &kase) {
-    REAL absxi = 0.0;
-    REAL altsgn = 0.0;
-    REAL estold = 0.0;
-    INTEGER i = 0;
-    INTEGER iter = 0;
-    INTEGER j = 0;
-    INTEGER jlast = 0;
-    INTEGER jump = 0;
-    REAL safmin = 0.0;
-    REAL temp = 0.0;
+    static REAL absxi = 0.0;
+    static REAL altsgn = 0.0;
+    static REAL estold = 0.0;
+    static INTEGER i = 0;
+    static INTEGER iter = 0;
+    static INTEGER j = 0;
+    static INTEGER jlast = 0;
+    static INTEGER jump = 0;
+    static REAL safmin = 0.0;
+    static REAL temp = 0.0;
     const REAL one = 1.0;
     const COMPLEX cone = COMPLEX(1.0, 0.0);
     const COMPLEX czero = COMPLEX(0.0, 0.0);
