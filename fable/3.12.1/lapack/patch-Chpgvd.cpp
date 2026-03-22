diff --git a/mplapack/reference/Chpgvd.cpp b/mplapack/reference/Chpgvd.cpp
index d0c3987b..868c9c5e 100644
--- a/mplapack/reference/Chpgvd.cpp
+++ b/mplapack/reference/Chpgvd.cpp
@@ -114,9 +114,9 @@ void Chpgvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     //
     Chpgst(itype, uplo, n, ap, bp, info);
     Chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
-    lwmin = castINTEGER(castREAL(max(castREAL(lwmin), work[1 - 1].real())));
-    lrwmin = castINTEGER(max(castREAL(lrwmin), castREAL(rwork[1 - 1])));
-    liwmin = castINTEGER(max(castREAL(liwmin), castREAL(iwork[1 - 1])));
+    lwmin = max(lwmin, castINTEGER(work[1 - 1].real()));
+    lrwmin = max(lrwmin, castINTEGER(rwork[1 - 1]));
+    liwmin = max(liwmin, iwork[1 - 1]);
     //
     INTEGER neig = 0;
     char trans;
