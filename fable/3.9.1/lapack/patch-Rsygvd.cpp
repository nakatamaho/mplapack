--- a/mplapack/reference/Rsygvd.cpp
+++ b/mplapack/reference/Rsygvd.cpp
@@ -109,8 +109,8 @@
     //
     Rsygst(itype, uplo, n, a, lda, b, ldb, info);
     Rsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
-    lopt = max(castREAL(lopt), castREAL(work[1 - 1]));
-    liopt = max(castREAL(liopt), castREAL(iwork[1 - 1]));
+    lopt = max(lopt, castINTEGER(work[1 - 1]));
+    liopt = max(liopt, iwork[1 - 1]);
     //
     char trans;
     const REAL one = 1.0;
