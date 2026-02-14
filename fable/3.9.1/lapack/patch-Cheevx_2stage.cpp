--- a/mplapack/reference/Cheevx_2stage.cpp
+++ b/mplapack/reference/Cheevx_2stage.cpp
@@ -273,6 +273,9 @@
     indisp = indibl + n;
     indiwk = indisp + n;
     Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &rwork[indd - 1], &rwork[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &rwork[indrwk - 1], &iwork[indiwk - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO from Rstebz; IBLOCK may be invalid
+    }
     //
     if (wantz) {
         Cstein(n, &rwork[indd - 1], &rwork[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &rwork[indrwk - 1], &iwork[indiwk - 1], ifail, info);
