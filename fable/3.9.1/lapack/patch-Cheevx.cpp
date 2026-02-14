https://claude.ai/share/5fa05feb-93e6-4a74-94e0-0a7b6aa0c6ab

see also MIGRATION_2.1.md
diff --git a/mplapack/reference/Cheevx.cpp b/mplapack/reference/Cheevx.cpp
index c72a90d1..f531f875 100644
--- a/mplapack/reference/Cheevx.cpp
+++ b/mplapack/reference/Cheevx.cpp
@@ -269,6 +269,21 @@ void Cheevx(const char *jobz, const char *range, const char *uplo, INTEGER const
     Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &rwork[indd - 1], &rwork[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &rwork[indrwk - 1], &iwork[indiwk - 1], info);
     //
     if (wantz) {
+        //
+        // When Rstebz reports non-convergence (info & 1), unconverged
+        // eigenvalues are marked with negative iblock values (iblock[j] = -jb).
+        // Cstein requires iblock to be monotonically non-decreasing and
+        // rejects negative entries with info = -6.  Skip Cstein in this case
+        // and let the Rstebz non-convergence info propagate to the caller.
+        //
+        if (info != 0) {
+            for (i = 1; i <= m; i = i + 1) {
+                if (iwork[(indibl + i - 1) - 1] < 0) {
+                    goto statement_40;
+                }
+            }
+        }
+        //
         Cstein(n, &rwork[indd - 1], &rwork[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &rwork[indrwk - 1], &iwork[indiwk - 1], ifail, info);
         //
         // Apply unitary matrix used in reduction to tridiagonal
