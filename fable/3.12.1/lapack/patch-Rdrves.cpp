--- Rdrves.cpp~	2026-01-22 14:21:54.586345859 +0900
+++ Rdrves.cpp	2026-01-22 14:38:42.624890439 +0900
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+
 void Rdrves(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER (&iseed)[4], REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *h, REAL *ht, REAL *wr, REAL *wi, REAL *wrt, REAL *wit, REAL *vs, INTEGER const ldvs, REAL *result, REAL *work, INTEGER const nwork, INTEGER *iwork, bool *bwork, INTEGER &info) {
     common cmn;
     common_write write(cmn);
