--- Rlctsx.cpp_	2026-01-22 13:39:22.837294653 +0900
+++ Rlctsx.cpp	2026-01-22 13:39:32.389459220 +0900
@@ -36,13 +36,12 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-#include <fem.hpp> // Fortran EMulation library of fable module
-using namespace fem::major_types;
-using fem::common;
-
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_mn.h>
+#include <mplapack_debug.h>
+
 bool Rlctsx(REAL const /* ar */, REAL const /* ai */, REAL const /* beta */) {
     bool return_value = false;
     //
