diff --git a/mplapack/reference/Rhgeqz.cpp b/mplapack/reference/Rhgeqz.cpp
index 666087997..b95f5a131 100644
--- a/mplapack/reference/Rhgeqz.cpp
+++ b/mplapack/reference/Rhgeqz.cpp
@@ -35,6 +35,8 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
+#include <mplapack_arithmetic_params_double.h>
 
 void Rhgeqz(const char *job, const char *compq, const char *compz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, REAL *h, INTEGER const ldh, REAL *t, INTEGER const ldt, REAL *alphar, REAL *alphai, REAL *beta, REAL *q, INTEGER const ldq, REAL *z, INTEGER const ldz, REAL *work, INTEGER const lwork, INTEGER &info) {
     bool ilschr = false;
@@ -63,8 +65,13 @@ void Rhgeqz(const char *job, const char *compq, const char *compz, INTEGER const
     INTEGER ilastm = 0;
     INTEGER iiter = 0;
     REAL eshift = 0.0;
+    const INTEGER maxit_base = 30;
     INTEGER maxit = 0;
+    const INTEGER maxit_cap = 1000;
     INTEGER jiter = 0;
+    REAL eps = 0.0;
+    REAL eps_ref = 0.0;
+    REAL scale_ratio = 0.0;
     bool ilazro = false;
     bool ilazr2 = false;
     REAL temp = 0.0;
@@ -240,7 +247,9 @@ void Rhgeqz(const char *job, const char *compq, const char *compz, INTEGER const
     in = ihi + 1 - ilo;
     safmin = Rlamch("S");
     safmax = one / safmin;
-    ulp = Rlamch("E") * Rlamch("B");
+    eps = Rlamch("E");
+    ulp = eps * Rlamch("B");
+    eps_ref = mplapack::get_arithmetic_params<double>().eps;
     anorm = Rlanhs("F", in, &h[(ilo - 1) + (ilo - 1) * ldh], ldh, work);
     bnorm = Rlanhs("F", in, &t[(ilo - 1) + (ilo - 1) * ldt], ldt, work);
     atol = max(safmin, ulp * anorm);
@@ -303,7 +312,16 @@ void Rhgeqz(const char *job, const char *compq, const char *compz, INTEGER const
     }
     iiter = 0;
     eshift = zero;
-    maxit = 30 * (ihi - ilo + 1);
+    // maxit is tuned for double precision. For high-precision backends,
+    // scale the QZ sweep limit with log(1/eps).
+    maxit = maxit_base * in;
+    scale_ratio = log(eps) / log(eps_ref);
+    if (scale_ratio >= one) {
+        maxit = iceil(castREAL(maxit_base * in) * scale_ratio);
+    }
+    if (maxit > maxit_cap) {
+        maxit = maxit_cap;
+    }
     //
     for (jiter = 1; jiter <= maxit; jiter = jiter + 1) {
         //
