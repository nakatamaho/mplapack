--- a/mplapack/reference/Chegvd.cpp
+++ b/mplapack/reference/Chegvd.cpp
@@ -117,9 +117,9 @@
     //
     Chegst(itype, uplo, n, a, lda, b, ldb, info);
     Cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
-    lopt = max(castREAL(lopt), work[1 - 1].real());
-    lropt = max(castREAL(lropt), castREAL(rwork[1 - 1]));
-    liopt = max(castREAL(liopt), castREAL(iwork[1 - 1]));
+    lopt = max(lopt, castINTEGER(work[1 - 1].real()));
+    lropt = max(lropt, castINTEGER(rwork[1 - 1]));
+    liopt = max(liopt, iwork[1 - 1]);
     //
     char trans;
     const COMPLEX cone = COMPLEX(1.0, 0.0);
