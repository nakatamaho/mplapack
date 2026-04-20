--- Rslect.cpp_	2026-01-22 13:36:11.902997305 +0900
+++ Rslect.cpp	2026-01-22 13:36:11.906997374 +0900
@@ -36,13 +36,12 @@
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
     bool return_value = false;
     //
