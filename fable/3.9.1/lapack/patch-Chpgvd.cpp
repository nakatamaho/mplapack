diff --git a/mplapack/reference/Chpgvd.cpp b/mplapack/reference/Chpgvd.cpp
index f68373fb..afa19de1 100644
--- a/mplapack/reference/Chpgvd.cpp
+++ b/mplapack/reference/Chpgvd.cpp
@@ -114,9 +114,9 @@ void Chpgvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     //
     Chpgst(itype, uplo, n, ap, bp, info);
     Chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
-    lwmin = max(castREAL(lwmin), work[0].real());
-    lrwmin = max(castREAL(lrwmin), castREAL(rwork[0]));
-    liwmin = max(castREAL(liwmin), castREAL(iwork[0]));
+    lwmin = max(lwmin, castINTEGER(work[0].real()));
+    lrwmin = max(lrwmin, castINTEGER(rwork[0]));
+    liwmin = max(liwmin, castINTEGER(iwork[0]));
     //
     INTEGER neig = 0;
     char trans;
