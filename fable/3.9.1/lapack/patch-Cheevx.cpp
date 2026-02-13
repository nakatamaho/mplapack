diff --git a/mplapack/reference/Cheevx.cpp b/mplapack/reference/Cheevx.cpp
index c72a90d1..cbf9d446 100644
--- a/mplapack/reference/Cheevx.cpp
+++ b/mplapack/reference/Cheevx.cpp
@@ -267,6 +267,9 @@ void Cheevx(const char *jobz, const char *range, const char *uplo, INTEGER const
     indisp = indibl + n;
     indiwk = indisp + n;
     Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &rwork[indd - 1], &rwork[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &rwork[indrwk - 1], &iwork[indiwk - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO (1/3) upwards
+    }
     //
     if (wantz) {
         Cstein(n, &rwork[indd - 1], &rwork[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &rwork[indrwk - 1], &iwork[indiwk - 1], ifail, info);
