--- a/mplapack/reference/Chpgvd.cpp
+++ b/mplapack/reference/Chpgvd.cpp
@@ -114,9 +114,9 @@
     //
     Chpgst(itype, uplo, n, ap, bp, info);
     Chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
-    lwmin = max(castREAL(lwmin), work[1 - 1].real());
-    lrwmin = max(castREAL(lrwmin), castREAL(rwork[1 - 1]));
-    liwmin = max(castREAL(liwmin), castREAL(iwork[1 - 1]));
+    lwmin = max(lwmin, castINTEGER(work[1 - 1].real()));
+    lrwmin = max(lrwmin, castINTEGER(rwork[1 - 1]));
+    liwmin = max(liwmin, iwork[1 - 1]);
     //
     INTEGER neig = 0;
     char trans;
