diff --git a/mplapack/test/eig/common/Cdrgsx.cpp b/mplapack/test/eig/common/Cdrgsx.cpp
index d39c39c9..1ae6c974 100644
--- Cdrgsx.cpp
+++ Cdrgsx.cpp
@@ -43,6 +43,8 @@ using fem::common;
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_mn.h>
+
 void Cdrgsx(INTEGER const nsize, INTEGER const ncmax, REAL const thresh, INTEGER const nin, INTEGER const nout, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *ai, COMPLEX *bi, COMPLEX *z, COMPLEX *q, COMPLEX *alpha, COMPLEX *beta, COMPLEX *c, INTEGER const ldc, REAL *s, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER const liwork, bool *bwork, INTEGER &info) {
     common cmn;
     common_read read(cmn);
@@ -254,7 +256,7 @@ void Cdrgsx(INTEGER const nsize, INTEGER const ncmax, REAL const thresh, INTEGER
                     //
                     for (j = 1; j <= mplusn; j = j + 1) {
                         ilabad = false;
-                        temp2 = (abs1(alpha[j - 1] - ai[(j - 1) + (j - 1) * lda]) / max(smlnum, abs1(alpha[j - 1]), abs1(ai[(j - 1) + (j - 1) * lda])) + abs1(beta[j - 1] - bi[(j - 1) + (j - 1) * lda]) / max(smlnum, abs1(beta[j - 1]), abs1(bi[(j - 1) + (j - 1) * lda]))) / ulp;
+                        temp2 = (cabs1(alpha[j - 1] - ai[(j - 1) + (j - 1) * lda]) / max(smlnum, cabs1(alpha[j - 1]), cabs1(ai[(j - 1) + (j - 1) * lda])) + cabs1(beta[j - 1] - bi[(j - 1) + (j - 1) * lda]) / max(smlnum, cabs1(beta[j - 1]), cabs1(bi[(j - 1) + (j - 1) * lda]))) / ulp;
                         if (j < mplusn) {
                             if (ai[((j + 1) - 1) + (j - 1) * lda] != zero) {
                                 ilabad = true;
@@ -472,7 +474,7 @@ statement_80:
     //
     for (j = 1; j <= mplusn; j = j + 1) {
         ilabad = false;
-        temp2 = (abs1(alpha[j - 1] - ai[(j - 1) + (j - 1) * lda]) / max(smlnum, abs1(alpha[j - 1]), abs1(ai[(j - 1) + (j - 1) * lda])) + abs1(beta[j - 1] - bi[(j - 1) + (j - 1) * lda]) / max(smlnum, abs1(beta[j - 1]), abs1(bi[(j - 1) + (j - 1) * lda]))) / ulp;
+        temp2 = (cabs1(alpha[j - 1] - ai[(j - 1) + (j - 1) * lda]) / max(smlnum, cabs1(alpha[j - 1]), cabs1(ai[(j - 1) + (j - 1) * lda])) + cabs1(beta[j - 1] - bi[(j - 1) + (j - 1) * lda]) / max(smlnum, cabs1(beta[j - 1]), cabs1(bi[(j - 1) + (j - 1) * lda]))) / ulp;
         if (j < mplusn) {
             if (ai[((j + 1) - 1) + (j - 1) * lda] != zero) {
                 ilabad = true;
