--- Cget24.cpp~	2026-01-22 14:21:56.444365288 +0900
+++ Cget24.cpp	2026-01-22 14:33:25.047309334 +0900
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+
 void Cget24(bool const comp, INTEGER const jtype, REAL const thresh, INTEGER (&iseed)[4], INTEGER const nounit, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *h, COMPLEX *ht, COMPLEX *w, COMPLEX *wt, COMPLEX *wtmp, COMPLEX *vs, INTEGER const ldvs, COMPLEX *vs1, REAL const rcdein, REAL const rcdvin, INTEGER const nslct, INTEGER *islct, INTEGER const isrt, REAL *result, COMPLEX *work, INTEGER const lwork, REAL *rwork, bool *bwork, INTEGER &info) {
     common cmn;
     common_write write(cmn);
