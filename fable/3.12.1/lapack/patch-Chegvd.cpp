diff --git a/mplapack/reference/Chegvd.cpp b/mplapack/reference/Chegvd.cpp
index deb4bc88..2b453edf 100644
--- a/mplapack/reference/Chegvd.cpp
+++ b/mplapack/reference/Chegvd.cpp
@@ -117,9 +117,9 @@ void Chegvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     //
     Chegst(itype, uplo, n, a, lda, b, ldb, info);
     Cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
-    lopt = castINTEGER(castREAL(max(castREAL(lopt), work[1 - 1].real())));
-    lropt = castINTEGER(max(castREAL(lropt), castREAL(rwork[1 - 1])));
-    liopt = castINTEGER(max(castREAL(liopt), castREAL(iwork[1 - 1])));
+    lopt = max(lopt, castINTEGER(work[1 - 1].real()));
+    lropt = max(lropt, castINTEGER(rwork[1 - 1]));
+    liopt = max(liopt, iwork[1 - 1]);
     //
     char trans;
     const COMPLEX cone = COMPLEX(1.0, 0.0);
