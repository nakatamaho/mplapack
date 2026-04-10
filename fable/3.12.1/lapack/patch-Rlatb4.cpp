--- Rlatb4.cpp	2026-04-10 17:38:16.783649501 +0900
+++ /home/docker/Rlatb4.cpp	2026-04-10 17:37:00.109473350 +0900
@@ -44,28 +44,33 @@
 #include <mplapack_lin.h>
 
 void Rlatb4(fem::str_cref path, INTEGER const imat, INTEGER const m, INTEGER const n, fem::str_ref type, INTEGER &kl, INTEGER &ku, REAL &anorm, INTEGER &mode, REAL &cndnum, fem::str_ref dist) {
-    static bool first = true;
-    double &badc1 = sve.badc1;
-    double &badc2 = sve.badc2;
-    double &eps = sve.eps;
-    double &large = sve.large;
-    double &small = sve.small;
     //
     // Set some constants for use in the subroutine.
     //
     const REAL tenth = 0.1;
     const REAL one = 1.0;
     const REAL shrink = 0.25;
-    if (first) {
-        first = false;
-        eps = Rlamch("Precision");
-        badc2 = tenth / eps;
-        badc1 = sqrt(badc2);
-        small = Rlamch("Safe minimum");
-        large = one / small;
-        small = shrink * (small / eps);
-        large = one / small;
+    REAL eps = Rlamch("Precision");
+    REAL badc2 = tenth / eps;
+    REAL badc1 = sqrt(badc2);
+    REAL small = Rlamch("Safe minimum");
+    REAL large = one / small;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+    // Historical DLABAD-style range reduction for very wide exponent ranges.
+    //
+    // Current LAPACK DLABAD is a no-op, but its former logic reduced SMALL and
+    // LARGE by square roots when LOG10(LARGE) > 2000. That condition is false
+    // for ordinary IEEE binary32/binary64 and true for wide-range formats such
+    // as binary80/binary128. Apply the same policy here so that standard
+    // backends keep LAPACK's original test scaling while wide-range backends
+    // avoid excessively extreme near-underflow/near-overflow matrices.
+    if (log10(large) > 2000.0) {
+        small = sqrt(small);
+        large = sqrt(large);
     }
+#endif
+    small = shrink * (small / eps);
+    large = one / small;
     //
     fem::str<2> c2 = path(2, 3);
     //
