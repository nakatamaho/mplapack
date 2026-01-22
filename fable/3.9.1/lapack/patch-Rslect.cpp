--- Rslect.cpp_	2026-01-22 08:18:09.082456791 +0900
+++ Rslect.cpp	2026-01-22 08:18:09.088456907 +0900
@@ -36,30 +36,25 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-#include <fem.hpp> // Fortran EMulation library of fable module
-using namespace fem::major_types;
-using fem::common;
-
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+#include <mplapack_debug.h>
+
 bool Rslect(REAL const zr, REAL const zi) {
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
         return_value = (zr < zero);
     } else {
         rmin = Rlapy2(zr - selwr[1 - 1], zi - selwi[1 - 1]);
         return_value = selval[1 - 1];
-        for (i = 2; i <= cmn.seldim; i = i + 1) {
+        for (i = 2; i <= seldim; i = i + 1) {
             x = Rlapy2(zr - selwr[i - 1], zi - selwi[i - 1]);
             if (x <= rmin) {
                 rmin = x;
