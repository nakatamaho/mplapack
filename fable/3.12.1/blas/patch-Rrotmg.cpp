--- Rrotmg.cpp_	2026-03-13 10:17:44.590377477 +0900
+++ Rrotmg.cpp	2026-03-13 10:17:48.618486134 +0900
@@ -36,12 +36,12 @@
 #include <mpblas.h>
 
 void Rrotmg(REAL &dd1, REAL &dd2, REAL &dx1, REAL const dy1, REAL *dparam) {
-    static REAL gam = 4096.0;
-    static REAL gamsq = 16777216.0;
-    static REAL one = 1.0;
-    static REAL rgamsq = 0.000000059604645;
-    static REAL two = 2.0;
-    static REAL zero = 0.0;
+    REAL gam = 4096.0;
+    REAL gamsq = 16777216.0;
+    REAL one = 1.0;
+    REAL rgamsq = 0x1p-24;
+    REAL two = 2.0;
+    REAL zero = 0.0;
     REAL dflag = 0.0;
     REAL dh11 = 0.0;
     REAL dh12 = 0.0;
