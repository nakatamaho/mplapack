diff --git a/mplapack/test/eig/common/Rdrgsx.cpp b/mplapack/test/eig/common/Rdrgsx.cpp
index 0a6d6d88..d5d9228e 100644
--- Rdrgsx.cpp
+++ Rdrgsx.cpp
@@ -43,6 +43,8 @@ using fem::common;
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_mn.h>
+
 void Rdrgsx(INTEGER const nsize, INTEGER const ncmax, REAL const thresh, INTEGER const nin, INTEGER const nout, REAL *a, INTEGER const lda, REAL *b, REAL *ai, REAL *bi, REAL *z, REAL *q, REAL *alphar, REAL *alphai, REAL *beta, REAL *c, INTEGER const ldc, REAL *s, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, bool *bwork, INTEGER &info) {
     common cmn;
     common_read read(cmn);
