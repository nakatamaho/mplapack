--- a/mplapack/reference/Rspgvd.cpp
+++ b/mplapack/reference/Rspgvd.cpp
@@ -106,8 +106,8 @@
     //
     Rspgst(itype, uplo, n, ap, bp, info);
     Rspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
-    lwmin = max(castREAL(lwmin), castREAL(work[1 - 1]));
-    liwmin = max(castREAL(liwmin), castREAL(iwork[1 - 1]));
+    lwmin = max(lwmin, castINTEGER(work[1 - 1]));
+    liwmin = max(liwmin, iwork[1 - 1]);
     //
     INTEGER neig = 0;
     char trans;
