diff --git a/mplapack/reference/Rsbgvx.cpp b/mplapack/reference/Rsbgvx.cpp
index 097cb646..feccf91f 100644
--- a/mplapack/reference/Rsbgvx.cpp
+++ b/mplapack/reference/Rsbgvx.cpp
@@ -188,6 +188,9 @@ void Rsbgvx(const char *jobz, const char *range, const char *uplo, INTEGER const
     indisp = indibl + n;
     indiwo = indisp + n;
     Rstebz(range, &order, n, vl, vu, il, iu, abstol, &work[indd - 1], &work[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &work[indwrk - 1], &iwork[indiwo - 1], info);
+    if (info != 0) {
+        return;  // propagate INFO (1/3) upwards
+    }
     //
     if (wantz) {
         Rstein(n, &work[indd - 1], &work[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &work[indwrk - 1], &iwork[indiwo - 1], ifail, info);
