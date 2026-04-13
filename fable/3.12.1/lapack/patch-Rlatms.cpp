--- Rlatms.cpp_	2026-01-17 18:56:26.497113566 +0900
+++ Rlatms.cpp	2026-01-17 18:56:34.394352153 +0900
@@ -275,7 +275,8 @@
     INTEGER jku = 0;
     INTEGER jr = 0;
     REAL extra = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL two = 2.0;
+    const REAL twopi = two * pi(zero);
     REAL angle = 0.0;
     REAL c = 0.0;
     REAL s = 0.0;

--- Rlatms.cpp
+++ Rlatms.cpp
@@ -329,7 +329,11 @@ void Rlatms(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&isee
                         extra = zero;
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        s = sqrt(one - c * c);
+#else
                         s = sin(angle);
+#endif
                         icol = max((INTEGER)1, jr - jkl);
                         if (jr < m) {
                             il = min(n, jr + jku) + 1 - icol;
@@ -371,7 +375,11 @@ void Rlatms(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&isee
                         extra = zero;
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        s = sqrt(one - c * c);
+#else
                         s = sin(angle);
+#endif
                         irow = max((INTEGER)1, jc - jku);
                         if (jc < n) {
                             il = min(m, jc + jkl) + 1 - irow;
@@ -421,7 +429,11 @@ void Rlatms(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&isee
                         extra = zero;
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        s = sqrt(one - c * c);
+#else
                         s = sin(angle);
+#endif
                         irow = max((INTEGER)1, jc - jku + 1);
                         if (jc > 0) {
                             il = min(m, jc + jkl + 1) + 1 - irow;
@@ -465,7 +477,11 @@ void Rlatms(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&isee
                         extra = zero;
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        s = sqrt(one - c * c);
+#else
                         s = sin(angle);
+#endif
                         icol = max((INTEGER)1, jr - jkl + 1);
                         if (jr > 0) {
                             il = min(n, jr + jku + 1) + 1 - icol;
@@ -524,7 +540,11 @@ void Rlatms(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&isee
                         temp = a[((jc - iskew * (jc + 1) + ioffg) - 1) + ((jc + 1) - 1) * lda];
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        s = sqrt(one - c * c);
+#else
                         s = sin(angle);
+#endif
                         Rlarot(false, jc > k, true, il, c, s, &a[((irow - iskew * jc + ioffg) - 1) + (jc - 1) * lda], ilda, extra, temp);
                         Rlarot(true, true, false, min(k, n - jc) + 1, c, s, &a[(((1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda], ilda, temp, dummy);
                         //
@@ -588,7 +608,11 @@ void Rlatms(INTEGER const m, INTEGER const n, fem::str_cref dist, INTEGER (&isee
                         temp = a[((1 + (1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda];
                         angle = twopi * Rlarnd(1, iseed);
                         c = cos(angle);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                        s = -sqrt(one - c * c);
+#else
                         s = -sin(angle);
+#endif
                         Rlarot(false, true, n - jc > k, il, c, s, &a[(((1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda], ilda, temp, extra);
                         icol = max((INTEGER)1, jc - k + 1);
                         Rlarot(true, false, true, jc + 2 - icol, c, s, &a[((jc - iskew * icol + ioffg) - 1) + (icol - 1) * lda], ilda, dummy, temp);
