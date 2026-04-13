diff --git a/mplapack/reference/Rbdsqr.cpp b/mplapack/reference/Rbdsqr.cpp
index c56dd7ce2..a241bb684 100644
--- a/mplapack/reference/Rbdsqr.cpp
+++ b/mplapack/reference/Rbdsqr.cpp
@@ -176,7 +176,12 @@ void Rbdsqr(const char *uplo, INTEGER const n, INTEGER const ncvt, INTEGER const
     // (By setting TOL to be negative, algorithm will compute
     // singular values to absolute accuracy ABS(TOL)*norm(input matrix))
     //
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+    tolmul = max(ten, min(hndrd, one / sqrt(sqrt(sqrt(eps)))));
+#else
     tolmul = max(ten, min(hndrd, pow(eps, meigth)));
+#endif
+
     tol = tolmul * eps;
     //
     // Compute approximate maximum, minimum singular values
