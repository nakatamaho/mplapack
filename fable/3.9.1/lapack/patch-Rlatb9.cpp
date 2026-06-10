--- Rlatb9.cpp_	2026-01-22 11:11:41.772820328 +0900
+++ Rlatb9.cpp	2026-01-22 11:12:02.972245459 +0900
@@ -44,33 +44,25 @@
 #include <mplapack_eig.h>
 
 void Rlatb9(fem::str_cref path, INTEGER const imat, INTEGER const m, INTEGER const p, INTEGER const n, fem::str_ref type, INTEGER &kla, INTEGER &kua, INTEGER &klb, INTEGER &kub, REAL &anorm, REAL &bnorm, INTEGER &modea, INTEGER &modeb, REAL &cndnma, REAL &cndnmb, fem::str_ref dista, fem::str_ref distb) {
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
-        //
-        // If it looks like we're on a Cray, take the square root of
-        // SMALL and LARGE to avoid overflow and underflow problems.
-        //
-        Rlabad(small, large);
-        small = shrink * (small / eps);
-        large = one / small;
-    }
+    REAL eps = Rlamch("Precision");
+    REAL badc2 = tenth / eps;
+    REAL badc1 = sqrt(badc2);
+    REAL small = Rlamch("Safe minimum");
+    REAL large = one / small;
+    //
+    // If it looks like we're on a Cray, take the square root of
+    // SMALL and LARGE to avoid overflow and underflow problems.
+    //
+    Rlabad(small, large);
+    small = shrink * (small / eps);
+    large = one / small;
+
     //
     // Set some parameters we don't plan to change.
     //
