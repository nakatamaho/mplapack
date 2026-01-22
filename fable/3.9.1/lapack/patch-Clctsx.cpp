--- Clctsx.cpp_	2026-01-22 13:36:11.874996820 +0900
+++ Clctsx.cpp	2026-01-22 13:36:11.881996941 +0900
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
 bool Clctsx(COMPLEX const /* alpha */, COMPLEX const /* beta */) {
     bool return_value = false;
     //
