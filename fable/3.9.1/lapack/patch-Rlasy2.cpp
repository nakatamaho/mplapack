--- a/mplapack/reference/Rlasy2.cpp
+++ b/mplapack/reference/Rlasy2.cpp
@@ -37,7 +37,6 @@
 #include <mplapack.h>
 
 void Rlasy2(bool const ltranl, bool const ltranr, INTEGER const isgn, INTEGER const n1, INTEGER const n2, REAL *tl, INTEGER const ldtl, REAL *tr, INTEGER const ldtr, REAL *b, INTEGER const ldb, REAL &scale, REAL *x, INTEGER const ldx, REAL &xnorm, INTEGER &info) {
-    common cmn;
     static INTEGER locu12[4] = {3, 4, 1, 2};
     static INTEGER locl21[4] = {2, 1, 4, 3};
     static INTEGER locu22[4] = {4, 3, 2, 1};
