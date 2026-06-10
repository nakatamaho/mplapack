--- Ccsdts.cpp_	2026-02-14 09:19:43.235825348 +0900
+++ Ccsdts.cpp	2026-02-14 09:19:48.811980317 +0900
@@ -184,8 +184,9 @@
     // Check sorting
     //
     const REAL realzero = 0.0;
+    const REAL two = 2.0;
     result[9 - 1] = realzero;
-    const REAL piover2 = 1.5707963267948966192313216916397514421;
+    const REAL piover2 = pi(realzero) / two;
     for (i = 1; i <= r; i = i + 1) {
         if (theta[i - 1] < realzero || theta[i - 1] > piover2) {
             result[9 - 1] = ulpinv;
