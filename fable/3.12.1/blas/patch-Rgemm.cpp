--- Rgemm.cpp.orig
+++ Rgemm.cpp
@@ -114,6 +114,9 @@
     //
     INTEGER l = 0;
     REAL temp = 0.0;
+#if defined MPLAPACK_BUILD_WITH_MPFR
+    temp.set_prec(mplapack_mpfr_rgemm_operation_precision(nota, notb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc));
+#endif
     if (notb) {
         if (nota) {
             //
             // Form  C := alpha*A*B + beta*C.
