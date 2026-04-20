--- Cdrgvx.cpp	2026-04-08 14:06:15
+++ Cdrgvx.cpp	2026-04-09 09:44:34
@@ -43,6 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+__attribute__((optimize("O1")))
 void Cdrgvx(INTEGER const nsize, REAL const thresh, INTEGER const nin, INTEGER const nout, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *ai, COMPLEX *bi, COMPLEX *alpha, COMPLEX *beta, COMPLEX *vl, COMPLEX *vr, INTEGER &ilo, INTEGER &ihi, REAL *lscale, REAL *rscale, REAL *s, REAL *dtru, REAL *dif, REAL *diftru, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER const liwork, REAL *result, bool *bwork, INTEGER &info) {
     common cmn;
     common_read read(cmn);
