--- Rlatmt.cpp_
+++ Rlatmt.cpp
@@ -38,6 +38,24 @@

 #include <mplapack_matgen.h>

+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+static void Rlatmt_random_unit_circle(INTEGER (&iseed)[4], REAL &c, REAL &s) {
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
 void Rlatmt(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&iseed)[4], fem::str_cref sym, REAL *d, INTEGER const mode, REAL const cond, REAL const dmax, INTEGER const rank, INTEGER const kl, INTEGER const ku, fem::str_cref pack, REAL *a, INTEGER const lda, REAL *work, INTEGER &info) {
     //
     // 1)      Decode and Test the input parameters.
@@ -275,7 +293,8 @@
     INTEGER jku = 0;
     INTEGER jr = 0;
     REAL extra = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL two = 2.0;
+    const REAL twopi = two * pi(zero);
     REAL angle = 0.0;
     REAL c = 0.0;
     REAL s = 0.0;
@@ -326,9 +345,13 @@
                     //
                     for (jr = 1; jr <= min(m + jku, n) + jkl - 1; jr = jr + 1) {
                         extra = zero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Rlatmt_random_unit_circle(iseed, c, s);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
                         s = sin(angle);
+#endif
                         icol = max((INTEGER)1, jr - jkl);
                         if (jr < m) {
                             il = min(n, jr + jku) + 1 - icol;
@@ -368,9 +391,13 @@
                     //
                     for (jc = 1; jc <= min(n + jkl, m) + jku - 1; jc = jc + 1) {
                         extra = zero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Rlatmt_random_unit_circle(iseed, c, s);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
                         s = sin(angle);
+#endif
                         irow = max((INTEGER)1, jc - jku);
                         if (jc < n) {
                             il = min(m, jc + jkl) + 1 - irow;
@@ -418,9 +445,13 @@
                     iendch = min(m, n + jkl) - 1;
                     for (jc = min(m + jku, n) - 1; jc >= 1 - jkl; jc = jc - 1) {
                         extra = zero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Rlatmt_random_unit_circle(iseed, c, s);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
                         s = sin(angle);
+#endif
                         irow = max((INTEGER)1, jc - jku + 1);
                         if (jc > 0) {
                             il = min(m, jc + jkl + 1) + 1 - irow;
@@ -462,9 +493,13 @@
                     iendch = min(n, m + jku) - 1;
                     for (jr = min(n + jkl, m) - 1; jr >= 1 - jku; jr = jr - 1) {
                         extra = zero;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Rlatmt_random_unit_circle(iseed, c, s);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
                         s = sin(angle);
+#endif
                         icol = max((INTEGER)1, jr - jkl + 1);
                         if (jr > 0) {
                             il = min(n, jr + jku + 1) + 1 - icol;
@@ -521,9 +556,13 @@
                         il = min(jc + 1, k + 2);
                         extra = zero;
                         temp = a[((jc - iskew * (jc + 1) + ioffg) - 1) + ((jc + 1) - 1) * lda];
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Rlatmt_random_unit_circle(iseed, c, s);
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
                         s = sin(angle);
+#endif
                         Rlarot(false, jc > k, true, il, c, s, &a[((irow - iskew * jc + ioffg) - 1) + (jc - 1) * lda], ilda, extra, temp);
                         Rlarot(true, true, false, min(k, n - jc) + 1, c, s, &a[(((1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda], ilda, temp, dummy);
                         //
@@ -585,9 +624,14 @@
                         il = min(n + 1 - jc, k + 2);
                         extra = zero;
                         temp = a[((1 + (1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda];
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        Rlatmt_random_unit_circle(iseed, c, s);
+                        s = -s;
+#else
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
                         s = -sin(angle);
+#endif
                         Rlarot(false, true, n - jc > k, il, c, s, &a[(((1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda], ilda, temp, extra);
                         icol = max((INTEGER)1, jc - k + 1);
                         Rlarot(true, false, true, jc + 2 - icol, c, s, &a[((jc - iskew * icol + ioffg) - 1) + (icol - 1) * lda], ilda, dummy, temp);
