diff --git a/mplapack/reference/Chegvd.cpp b/mplapack/reference/Chegvd.cpp
index daead295..81a24bac 100644
--- a/mplapack/reference/Chegvd.cpp
+++ b/mplapack/reference/Chegvd.cpp
@@ -117,9 +117,9 @@ void Chegvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     //
     Chegst(itype, uplo, n, a, lda, b, ldb, info);
     Cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
-    lopt = max(castREAL(lopt), work[0].real());
-    lropt = max(castREAL(lropt), castREAL(rwork[0]));
-    liopt = max(castREAL(liopt), castREAL(iwork[0]));
+    lopt = max(lopt, castINTEGER(work[0].real()));
+    lropt = max(lropt, castINTEGER(rwork[0]));
+    liopt = max(liopt, iwork[0]);
     //
     char trans;
     const COMPLEX cone = COMPLEX(1.0, 0.0);
