diff --git a/mplapack/reference/Chpgvd.cpp b/mplapack/reference/Chpgvd.cpp
index afa19de1..81cc37fe 100644
--- a/mplapack/reference/Chpgvd.cpp
+++ b/mplapack/reference/Chpgvd.cpp
@@ -116,7 +116,7 @@ void Chpgvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     Chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
     lwmin = max(lwmin, castINTEGER(work[0].real()));
     lrwmin = max(lrwmin, castINTEGER(rwork[0]));
-    liwmin = max(liwmin, castINTEGER(iwork[0]));
+    liwmin = max(liwmin, iwork[0]);
     //
     INTEGER neig = 0;
     char trans;
