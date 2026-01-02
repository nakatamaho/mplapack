--- a/mplapack/reference/Clatbs.cpp
+++ b/mplapack/reference/Clatbs.cpp
@@ -36,6 +36,8 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL cabs2(COMPLEX zdum) { return abs(zdum.real() / 2.0) + abs(zdum.imag() / 2.0); }
+
 void Clatbs(const char *uplo, const char *trans, const char *diag, const char *normin, INTEGER const n, INTEGER const kd, COMPLEX *ab, INTEGER const ldab, COMPLEX *x, REAL &scale, REAL *cnorm, INTEGER &info) {
     COMPLEX zdum = 0.0;
     bool upper = false;
