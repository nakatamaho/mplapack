diff --git a/mplapack/reference/Rstebz.cpp b/mplapack/reference/Rstebz.cpp
index c16fca96..47f4ffeb 100644
--- a/mplapack/reference/Rstebz.cpp
+++ b/mplapack/reference/Rstebz.cpp
@@ -223,10 +223,6 @@ void Rstebz(const char *range, const char *order, INTEGER const n, REAL const vl
         // Compute Iteration parameters
         //
         itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
-#if defined ___MPLAPACK_BUILD_WITH_MPFR___
-        if (itmax >= 10240)
-            itmax = 10240; // XXX itmax can be too large for MPFR (=10^8)
-#endif
         if (abstol <= zero) {
             atoli = ulp * tnorm;
         } else {
@@ -382,10 +378,6 @@ void Rstebz(const char *range, const char *order, INTEGER const n, REAL const vl
             // Compute Eigenvalues
             //
             itmax = castINTEGER((log(gu - gl + pivmin) - log(pivmin)) / log(two)) + 2;
-#if defined ___MPLAPACK_BUILD_WITH_MPFR___
-            if (itmax >= 10240)
-                itmax = 10240; // XXX itmax can be too large for MPFR (=10^8)
-#endif
             Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &work[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], iout, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
             //
             // Copy Eigenvalues Into W and IBLOCK
