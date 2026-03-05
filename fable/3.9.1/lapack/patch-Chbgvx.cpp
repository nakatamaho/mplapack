--- a/mplapack/reference/Chbgvx.cpp
+++ b/mplapack/reference/Chbgvx.cpp
@@ -192,6 +192,9 @@
     indisp = indibl + n;
     indiwk = indisp + n;
     Rstebz(range, &order, n, vl, vu, il, iu, abstol, &rwork[indd - 1], &rwork[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &rwork[indrwk - 1], &iwork[indiwk - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO from Rstebz; IBLOCK may be invalid
+    }
     //
     if (wantz) {
         Cstein(n, &rwork[indd - 1], &rwork[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &rwork[indrwk - 1], &iwork[indiwk - 1], ifail, info);
