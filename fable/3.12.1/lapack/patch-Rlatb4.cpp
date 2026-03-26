--- Rlatb4.cpp_	2026-03-23 17:33:50.421604561 +0900
+++ Rlatb4.cpp	2026-03-23 17:33:57.295751242 +0900
@@ -44,28 +44,19 @@
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
-    }
+    REAL eps = Rlamch("Precision");
+    REAL badc2 = tenth / eps;
+    REAL badc1 = sqrt(badc2);
+    REAL small = Rlamch("Safe minimum");
+    REAL large = one / small;
+    small = shrink * (small / eps);
+    large = one / small;
     //
     fem::str<2> c2 = path(2, 3);
     //
--- Rlatb4.cpp	2026-03-26 18:16:03.465348200 +0900
+++ Rlatb4.cpp	2026-03-26 17:59:18.957020593 +0900
@@ -55,6 +55,18 @@
     REAL badc1 = sqrt(badc2);
     REAL small = Rlamch("Safe minimum");
     REAL large = one / small;
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
+    }
     small = shrink * (small / eps);
     large = one / small;
     //
