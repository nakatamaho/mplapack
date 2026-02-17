--- Rget24.cpp~	2026-01-22 14:21:54.724347300 +0900
+++ Rget24.cpp	2026-01-22 14:37:07.985193701 +0900
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+
 void Rget24(bool const comp, INTEGER const jtype, REAL const thresh, INTEGER (&iseed)[4], INTEGER const nounit, INTEGER const n, REAL *a, INTEGER const lda, REAL *h, REAL *ht, REAL *wr, REAL *wi, REAL *wrt, REAL *wit, REAL *wrtmp, REAL *witmp, REAL *vs, INTEGER const ldvs, REAL *vs1, REAL const rcdein, REAL const rcdvin, INTEGER const nslct, INTEGER *islct, REAL *result, REAL *work, INTEGER const lwork, INTEGER *iwork, bool *bwork, INTEGER &info) {
     common cmn;
     common_write write(cmn);
