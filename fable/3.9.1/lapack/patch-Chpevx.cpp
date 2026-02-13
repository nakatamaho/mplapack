diff --git a/mplapack/reference/Chpevx.cpp b/mplapack/reference/Chpevx.cpp
index b7a84b71..a8bab8e5 100644
--- a/mplapack/reference/Chpevx.cpp
+++ b/mplapack/reference/Chpevx.cpp
@@ -231,6 +231,9 @@ void Chpevx(const char *jobz, const char *range, const char *uplo, INTEGER const
     indisp = indibl + n;
     indiwk = indisp + n;
     Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &rwork[indd - 1], &rwork[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &rwork[indrwk - 1], &iwork[indiwk - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO (1/3) upwards
+    }
     //
     if (wantz) {
         Cstein(n, &rwork[indd - 1], &rwork[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &rwork[indrwk - 1], &iwork[indiwk - 1], ifail, info);
