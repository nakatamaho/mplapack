diff --git a/mplapack/reference/Rsyevx.cpp b/mplapack/reference/Rsyevx.cpp
index 4a58a75c..a786762a 100644
--- a/mplapack/reference/Rsyevx.cpp
+++ b/mplapack/reference/Rsyevx.cpp
@@ -266,6 +266,9 @@ void Rsyevx(const char *jobz, const char *range, const char *uplo, INTEGER const
     indisp = indibl + n;
     indiwo = indisp + n;
     Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &work[indd - 1], &work[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &work[indwrk - 1], &iwork[indiwo - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO (1/3) upwards
+    }
     //
     if (wantz) {
         Rstein(n, &work[indd - 1], &work[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &work[indwrk - 1], &iwork[indiwo - 1], ifail, info);
