--- Cslect.cpp_	2026-01-22 08:18:09.053456228 +0900
+++ Cslect.cpp	2026-01-22 08:18:09.059456345 +0900
@@ -43,23 +43,22 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+#include <mplapack_debug.h>
+
 bool Cslect(COMPLEX const z) {
-    common cmn;
     bool return_value = false;
-    arr_cref<bool> selval(cmn.selval, dimension(20));
-    arr_cref<double> selwr(cmn.selwr, dimension(20));
-    arr_cref<double> selwi(cmn.selwi, dimension(20));
     //
     const REAL zero = 0.0;
     REAL rmin = 0.0;
     INTEGER i = 0;
     REAL x = 0.0;
-    if (cmn.selopt == 0) {
+    if (selopt == 0) {
         return_value = (z.real() < zero);
     } else {
         rmin = abs(z - COMPLEX(selwr[1 - 1], selwi[1 - 1]));
         return_value = selval[1 - 1];
-        for (i = 2; i <= cmn.seldim; i = i + 1) {
+        for (i = 2; i <= seldim; i = i + 1) {
             x = abs(z - COMPLEX(selwr[i - 1], selwi[i - 1]));
             if (x <= rmin) {
                 rmin = x;
