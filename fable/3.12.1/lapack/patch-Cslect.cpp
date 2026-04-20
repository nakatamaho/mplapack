--- Cslect.cpp_	2026-01-22 13:36:11.889997079 +0900
+++ Cslect.cpp	2026-01-22 13:36:11.893997149 +0900
@@ -43,6 +43,9 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
+#include <mplapack_common_sslct.h>
+#include <mplapack_debug.h>
+
 bool Cslect(COMPLEX const z) {
     bool return_value = false;
     //
