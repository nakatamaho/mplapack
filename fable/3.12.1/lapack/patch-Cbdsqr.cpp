diff --git a/mplapack/reference/Cbdsqr.cpp b/mplapack/reference/Cbdsqr.cpp
index afc1779af..40be7bc43 100644
--- a/mplapack/reference/Cbdsqr.cpp
+++ b/mplapack/reference/Cbdsqr.cpp
@@ -36,2 +36,4 @@
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
+#include <mplapack_arithmetic_params_double.h>
@@ -59,7 +59,11 @@ void Cbdsqr(const char *uplo, INTEGER const n, INTEGER const ncvt, INTEGER const
     REAL smin = 0.0;
     REAL sminoa = 0.0;
     REAL mu = 0.0;
-    const INTEGER maxitr = 6;
+    const INTEGER maxitr_base = 6;
+    INTEGER maxitr = maxitr_base;
+    const INTEGER maxitr_cap = 1000;
+    REAL eps_ref = 0.0;
+    REAL scale_ratio = 0.0;
     REAL thresh = 0.0;
     INTEGER maxitdivn = 0;
     INTEGER iterdivn = 0;
@@ -148,6 +152,16 @@ void Cbdsqr(const char *uplo, INTEGER const n, INTEGER const ncvt, INTEGER const
     //
     eps = Rlamch("Epsilon");
     unfl = Rlamch("Safe minimum");
+    eps_ref = mplapack::get_arithmetic_params<double>().eps;
+    // maxitr is tuned for double precision. For high-precision backends,
+    // scale the zero-shift QR sweep limit with log(1/eps).
+    scale_ratio = log(eps) / log(eps_ref);
+    if (scale_ratio >= one) {
+        maxitr = iceil(castREAL(maxitr_base) * scale_ratio);
+    }
+    if (maxitr > maxitr_cap) {
+        maxitr = maxitr_cap;
+    }
     //
     // If matrix lower bidiagonal, rotate to be upper bidiagonal
     // by applying Givens rotations on the left
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
