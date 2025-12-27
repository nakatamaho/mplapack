diff --git a/mplapack/reference/Rspgvd.cpp b/mplapack/reference/Rspgvd.cpp
index 57ac111b..15a81f33 100644
--- a/mplapack/reference/Rspgvd.cpp
+++ b/mplapack/reference/Rspgvd.cpp
@@ -106,8 +106,8 @@ void Rspgvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER con
     //
     Rspgst(itype, uplo, n, ap, bp, info);
     Rspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
-    lwmin = max(castREAL(lwmin), castREAL(work[0]));
-    liwmin = max(castREAL(liwmin), castREAL(iwork[0]));
+    lwmin = max(lwmin, castINTEGER(work[0]));
+    liwmin = max(liwmin, iwork[0]);
     //
     INTEGER neig = 0;
     char trans;
