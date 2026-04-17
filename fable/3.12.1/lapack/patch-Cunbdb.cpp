--- a/mplapack/reference/Cunbdb.cpp
+++ b/mplapack/reference/Cunbdb.cpp
@@ -37,6 +37,10 @@
 #include <mplapack.h>
 
 void Cunbdb(const char *trans, const char *signs, INTEGER const m, INTEGER const p, INTEGER const q, COMPLEX *x11, INTEGER const ldx11, COMPLEX *x12, INTEGER const ldx12, COMPLEX *x21, INTEGER const ldx21, COMPLEX *x22, INTEGER const ldx22, REAL *theta, REAL *phi, COMPLEX *taup1, COMPLEX *taup2, COMPLEX *tauq1, COMPLEX *tauq2, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+    printf("MPLAPACK ERROR Cunbdb.cpp is not supported for GMP\n");
+    exit(1);
+#endif
     //
     // Test input arguments
     //
