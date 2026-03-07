--- a/mplapack/reference/Clatrs.cpp
+++ b/mplapack/reference/Clatrs.cpp
@@ -36,6 +36,8 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL cabs2(COMPLEX zdum) { return abs(zdum.real() / 2.0) + abs(zdum.imag() / 2.0); }
+
 void Clatrs(const char *uplo, const char *trans, const char *diag, const char *normin, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, REAL &scale, REAL *cnorm, INTEGER &info) {
     COMPLEX zdum = 0.0;
     bool upper = false;
