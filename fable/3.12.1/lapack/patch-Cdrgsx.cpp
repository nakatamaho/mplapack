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
