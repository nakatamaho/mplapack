--- a/mplapack/reference/Rlacon.cpp
+++ b/mplapack/reference/Rlacon.cpp
@@ -37,15 +37,14 @@
 #include <mplapack.h>
 
 void Rlacon(INTEGER const n, REAL *v, REAL *x, INTEGER *isgn, REAL &est, INTEGER &kase) {
-    common cmn;
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
