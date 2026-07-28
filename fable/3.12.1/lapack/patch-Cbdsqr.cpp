--- a/mplapack/reference/Cbdsqr.cpp
+++ b/mplapack/reference/Cbdsqr.cpp
@@ -35,6 +35,8 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
+#include <mplapack_arithmetic_params_double.h>
 
 void Cbdsqr(const char *uplo, INTEGER const n, INTEGER const ncvt, INTEGER const nru, INTEGER const ncc, REAL *d, REAL *e, COMPLEX *vt, INTEGER const ldvt, COMPLEX *u, INTEGER const ldu, COMPLEX *c, INTEGER const ldc, REAL *rwork, INTEGER &info) {
     bool lower = false;
@@ -59,7 +61,11 @@
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
@@ -148,6 +154,16 @@
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
@@ -176,7 +192,11 @@
     // (By setting TOL to be negative, algorithm will compute
     // singular values to absolute accuracy ABS(TOL)*norm(input matrix))
     //
+#if defined MPLAPACK_BUILD_WITH_GMP
+    tolmul = max(ten, min(hndrd, one / sqrt(sqrt(sqrt(eps)))));
+#else
     tolmul = max(ten, min(hndrd, pow(eps, meigth)));
+#endif
     tol = tolmul * eps;
     //
     // Compute approximate maximum, minimum singular values
