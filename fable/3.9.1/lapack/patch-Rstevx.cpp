--- a/mplapack/reference/Rstevx.cpp
+++ b/mplapack/reference/Rstevx.cpp
@@ -208,6 +208,9 @@
     indisp = indibl + n;
     indiwo = indisp + n;
     Rstebz(range, &order, n, vll, vuu, il, iu, abstol, d, e, m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &work[indwrk - 1], &iwork[indiwo - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO from Rstebz; IBLOCK may be invalid
+    }
     //
     if (wantz) {
         Rstein(n, d, e, m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &work[indwrk - 1], &iwork[indiwo - 1], ifail, info);
