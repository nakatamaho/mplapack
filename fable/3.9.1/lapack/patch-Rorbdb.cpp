diff --git a/mplapack/reference/Rorbdb.cpp b/mplapack/reference/Rorbdb.cpp
index aead7009..727cb2a8 100644
--- a/mplapack/reference/Rorbdb.cpp
+++ b/mplapack/reference/Rorbdb.cpp
@@ -37,6 +37,10 @@
 #include <mplapack.h>
 
 void Rorbdb(const char *trans, const char *signs, INTEGER const m, INTEGER const p, INTEGER const q, REAL *x11, INTEGER const ldx11, REAL *x12, INTEGER const ldx12, REAL *x21, INTEGER const ldx21, REAL *x22, INTEGER const ldx22, REAL *theta, REAL *phi, REAL *taup1, REAL *taup2, REAL *tauq1, REAL *tauq2, REAL *work, INTEGER const lwork, INTEGER &info) {
+#if defined MPLAPACK_BUILD_WITH_GMP
+    printf("MPLAPACK ERROR Rorbdb.cpp is not supported for GMP\n");
+    exit(1);
+#endif
     //
     // Test input arguments
     //
