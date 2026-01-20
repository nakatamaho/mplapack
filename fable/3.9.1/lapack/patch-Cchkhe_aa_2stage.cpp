--- Cchkhe_aa_2stage.cpp
+++ Cchkhe_aa_2stage.cpp
@@ -43,7 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void Cchkhe_aa_2stage(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const nmax, COMPLEX *a, COMPLEX *afac, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
+void Cchkhe_aa_2stage(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const nmax, COMPLEX *a, COMPLEX *afac, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, COMPLEX *rwork, INTEGER *iwork, INTEGER const nout) {
     common cmn;
     common_write write(cmn);
     static INTEGER iseedy[4] = {1988, 1989, 1990, 1991};
