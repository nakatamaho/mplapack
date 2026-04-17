--- a/mplapack/reference/Chgeqz.cpp
+++ b/mplapack/reference/Chgeqz.cpp
@@ -35,6 +35,8 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
+#include <mplapack_arithmetic_params_double.h>
 
 void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, COMPLEX *h, INTEGER const ldh, COMPLEX *t, INTEGER const ldt, COMPLEX *alpha, COMPLEX *beta, COMPLEX *q, INTEGER const ldq, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER &info) {
     COMPLEX x = 0.0;
@@ -65,8 +67,13 @@
     INTEGER ilastm = 0;
     INTEGER iiter = 0;
     COMPLEX eshift = 0.0;
+    const INTEGER maxit_base = 30;
     INTEGER maxit = 0;
+    const INTEGER maxit_cap = 1000;
     INTEGER jiter = 0;
+    REAL eps = 0.0;
+    REAL eps_ref = 0.0;
+    REAL scale_ratio = 0.0;
     bool ilazro = false;
     bool ilazr2 = false;
     INTEGER jch = 0;
@@ -191,7 +198,9 @@
     //
     in = ihi + 1 - ilo;
     safmin = Rlamch("S");
-    ulp = Rlamch("E") * Rlamch("B");
+    eps = Rlamch("E");
+    ulp = eps * Rlamch("B");
+    eps_ref = mplapack::get_arithmetic_params<double>().eps;
     anorm = Clanhs("F", in, &h[(ilo - 1) + (ilo - 1) * ldh], ldh, rwork);
     bnorm = Clanhs("F", in, &t[(ilo - 1) + (ilo - 1) * ldt], ldt, rwork);
     atol = max(safmin, ulp * anorm);
@@ -253,7 +262,16 @@
     }
     iiter = 0;
     eshift = czero;
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
