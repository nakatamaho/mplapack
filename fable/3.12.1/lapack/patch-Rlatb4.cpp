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
