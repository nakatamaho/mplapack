diff --git a/mplapack/reference/Cbdsqr.cpp b/mplapack/reference/Cbdsqr.cpp
index afc1779af..85542972a 100644
--- a/mplapack/reference/Cbdsqr.cpp
+++ b/mplapack/reference/Cbdsqr.cpp
@@ -176,7 +176,11 @@ void Cbdsqr(const char *uplo, INTEGER const n, INTEGER const ncvt, INTEGER const
     // (By setting TOL to be negative, algorithm will compute
     // singular values to absolute accuracy ABS(TOL)*norm(input matrix))
     //
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+    tolmul = max(ten, min(hndrd, one / sqrt(sqrt(sqrt(eps)))));
+#else
     tolmul = max(ten, min(hndrd, pow(eps, meigth)));
+#endif
     tol = tolmul * eps;
     //
     // Compute approximate maximum, minimum singular values
