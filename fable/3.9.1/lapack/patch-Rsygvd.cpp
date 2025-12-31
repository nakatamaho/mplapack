diff --git a/mplapack/reference/Rsygvd.cpp b/mplapack/reference/Rsygvd.cpp
index 3b1f6f00..b030c35d 100644
--- a/mplapack/reference/Rsygvd.cpp
+++ b/mplapack/reference/Rsygvd.cpp
@@ -109,8 +109,8 @@ void Rsygvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     //
     Rsygst(itype, uplo, n, a, lda, b, ldb, info);
     Rsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
-    lopt = max(castREAL(lopt), castREAL(work[0]));
-    liopt = max(castREAL(liopt), castREAL(iwork[0]));
+    lopt = max(lopt, castINTEGER(work[0]));
+    liopt = max(liopt, iwork[0]);
     //
     char trans;
     const REAL one = 1.0;
