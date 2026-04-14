--- Clatmt.cpp
+++ Clatmt.cpp
@@ -38,6 +38,24 @@

 #include <mplapack_matgen.h>

+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+static void Clatmt_random_unit_circle(INTEGER (&iseed)[4], REAL &c, REAL &s) {
+    const REAL zero = 0.0;
+    const REAL one = 1.0;
+    const REAL two = 2.0;
+    REAL u = 0.0;
+    REAL v = 0.0;
+    REAL r = 0.0;
+    do {
+        u = two * Rlarnd(1, iseed) - one;
+        v = two * Rlarnd(1, iseed) - one;
+        r = u * u + v * v;
+    } while (r <= zero || r > one);
+    c = (u * u - v * v) / r;
+    s = two * u * v / r;
+}
+#endif
+
 void Clatmt(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&iseed)[4], fem::str_cref sym, REAL *d, INTEGER const mode, REAL const cond, REAL const dmax, INTEGER const rank, INTEGER const kl, INTEGER const ku, fem::str_cref pack, COMPLEX *a, INTEGER const lda, COMPLEX *work, INTEGER &info) {
     //
     // 1)      Decode and Test the input parameters.
@@ -283,7 +301,8 @@
     INTEGER jku = 0;
     INTEGER jr = 0;
     COMPLEX extra = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL two = 2.0;
+    const REAL twopi = two * pi(zero);
     REAL angle = 0.0;
     COMPLEX c = 0.0;
     COMPLEX s = 0.0;
@@ -294,6 +313,10 @@
     INTEGER ic = 0;
     INTEGER jch = 0;
     REAL realc = 0.0;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+    REAL cr = 0.0;
+    REAL sr = 0.0;
+#endif
     INTEGER irow = 0;
     COMPLEX ztemp = 0.0;
     bool iltemp = false;
@@ -343,9 +366,15 @@
                     //
                     for (jr = 1; jr <= min(m + jku, n) + jkl - 1; jr = jr + 1) {
                         extra = czero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Clatmt_random_unit_circle(iseed, cr, sr);
+                        c = cr * Clarnd(5, iseed);
+                        s = sr * Clarnd(5, iseed);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle) * Clarnd(5, iseed);
                         s = sin(angle) * Clarnd(5, iseed);
+#endif
                         icol = max((INTEGER)1, jr - jkl);
                         if (jr < m) {
                             il = min(n, jr + jku) + 1 - icol;
@@ -392,9 +421,15 @@
                     //
                     for (jc = 1; jc <= min(n + jkl, m) + jku - 1; jc = jc + 1) {
                         extra = czero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Clatmt_random_unit_circle(iseed, cr, sr);
+                        c = cr * Clarnd(5, iseed);
+                        s = sr * Clarnd(5, iseed);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle) * Clarnd(5, iseed);
                         s = sin(angle) * Clarnd(5, iseed);
+#endif
                         irow = max((INTEGER)1, jc - jku);
                         if (jc < n) {
                             il = min(m, jc + jkl) + 1 - irow;
@@ -448,9 +483,15 @@
                     iendch = min(m, n + jkl) - 1;
                     for (jc = min(m + jku, n) - 1; jc >= 1 - jkl; jc = jc - 1) {
                         extra = czero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Clatmt_random_unit_circle(iseed, cr, sr);
+                        c = cr * Clarnd(5, iseed);
+                        s = sr * Clarnd(5, iseed);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle) * Clarnd(5, iseed);
                         s = sin(angle) * Clarnd(5, iseed);
+#endif
                         irow = max((INTEGER)1, jc - jku + 1);
                         if (jc > 0) {
                             il = min(m, jc + jkl + 1) + 1 - irow;
@@ -498,9 +539,15 @@
                     iendch = min(n, m + jku) - 1;
                     for (jr = min(n + jkl, m) - 1; jr >= 1 - jku; jr = jr - 1) {
                         extra = czero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Clatmt_random_unit_circle(iseed, cr, sr);
+                        c = cr * Clarnd(5, iseed);
+                        s = sr * Clarnd(5, iseed);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle) * Clarnd(5, iseed);
                         s = sin(angle) * Clarnd(5, iseed);
+#endif
                         icol = max((INTEGER)1, jr - jkl + 1);
                         if (jr > 0) {
                             il = min(n, jr + jku + 1) + 1 - icol;
@@ -568,9 +615,15 @@
                         il = min(jc + 1, k + 2);
                         extra = czero;
                         ztemp = a[((jc - iskew * (jc + 1) + ioffg) - 1) + ((jc + 1) - 1) * lda];
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Clatmt_random_unit_circle(iseed, cr, sr);
+                        c = cr * Clarnd(5, iseed);
+                        s = sr * Clarnd(5, iseed);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle) * Clarnd(5, iseed);
                         s = sin(angle) * Clarnd(5, iseed);
+#endif
                         if (csym) {
                             ct = c;
                             st = s;
@@ -660,9 +713,15 @@
                         il = min(n + 1 - jc, k + 2);
                         extra = czero;
                         ztemp = a[((1 + (1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda];
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Clatmt_random_unit_circle(iseed, cr, sr);
+                        c = cr * Clarnd(5, iseed);
+                        s = sr * Clarnd(5, iseed);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle) * Clarnd(5, iseed);
                         s = sin(angle) * Clarnd(5, iseed);
+#endif
                         if (csym) {
                             ct = c;
                             st = s;
