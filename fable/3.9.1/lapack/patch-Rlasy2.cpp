--- a/mplapack/reference/Rlasy2.cpp
+++ b/mplapack/reference/Rlasy2.cpp
@@ -37,11 +37,11 @@
 #include <mplapack.h>
 
 void Rlasy2(bool const ltranl, bool const ltranr, INTEGER const isgn, INTEGER const n1, INTEGER const n2, REAL *tl, INTEGER const ldtl, REAL *tr, INTEGER const ldtr, REAL *b, INTEGER const ldb, REAL &scale, REAL *x, INTEGER const ldx, REAL &xnorm, INTEGER &info) {
-    bool bswpiv[4];
-    INTEGER locl21[4];
-    INTEGER locu12[4];
-    INTEGER locu22[4];
-    bool xswpiv[4];
+    static INTEGER locu12[4] = {3, 4, 1, 2};
+    static INTEGER locl21[4] = {2, 1, 4, 3};
+    static INTEGER locu22[4] = {4, 3, 2, 1};
+    static bool xswpiv[4] = {false, false, true, true};
+    static bool bswpiv[4] = {false, true, false, true};
     REAL eps = 0.0;
     REAL smlnum = 0.0;
     REAL sgn = 0.0;
