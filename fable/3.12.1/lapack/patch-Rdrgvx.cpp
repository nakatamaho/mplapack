--- Rdrgvx.cpp	2026-04-08 13:57:12
+++ Rdrgvx.cpp	2026-04-08 13:53:04
@@ -43,6 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+__attribute__((optimize("O1")))
 void Rdrgvx(INTEGER const nsize, REAL const thresh, INTEGER const nin, INTEGER const nout, REAL *a, INTEGER const lda, REAL *b, REAL *ai, REAL *bi, REAL *alphar, REAL *alphai, REAL *beta, REAL *vl, REAL *vr, INTEGER &ilo, INTEGER &ihi, REAL *lscale, REAL *rscale, REAL *s, REAL *dtru, REAL *dif, REAL *diftru, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, REAL *result, bool *bwork, INTEGER &info) {
     common cmn;
     common_read read(cmn);
