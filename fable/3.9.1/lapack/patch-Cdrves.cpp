--- Cdrves.cpp~	2026-01-22 14:21:56.187362597 +0900
+++ Cdrves.cpp	2026-01-22 14:35:13.716180215 +0900
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+
 void Cdrves(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER (&iseed)[4], REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *h, COMPLEX *ht, COMPLEX *w, COMPLEX *wt, COMPLEX *vs, INTEGER const ldvs, REAL *result, COMPLEX *work, INTEGER const nwork, REAL *rwork, INTEGER *iwork, bool *bwork, INTEGER &info) {
     common cmn;
     common_write write(cmn);
