--- a/mplapack/reference/Rsbevx_2stage.cpp
+++ b/mplapack/reference/Rsbevx_2stage.cpp
@@ -273,6 +273,9 @@
     indisp = indibl + n;
     indiwo = indisp + n;
     Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &work[indd - 1], &work[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &work[indwrk - 1], &iwork[indiwo - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO from Rstebz; IBLOCK may be invalid
+    }
     //
     if (wantz) {
         Rstein(n, &work[indd - 1], &work[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &work[indwrk - 1], &iwork[indiwo - 1], ifail, info);
