--- Rlctsx.cpp_	2026-01-22 08:18:09.067456500 +0900
+++ Rlctsx.cpp	2026-01-22 08:18:09.072456597 +0900
@@ -36,23 +36,19 @@
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
     common cmn;
     bool return_value = false;
-    int &mplusn = cmn.mplusn;
-    int &i = cmn.i;
-    bool &fs = cmn.fs;
     //
     if (fs) {
         i++;
-        if (i <= cmn.m) {
+        if (i <= m) {
             return_value = false;
         } else {
             return_value = true;
@@ -63,7 +59,7 @@
         }
     } else {
         i++;
-        if (i <= cmn.n) {
+        if (i <= n) {
             return_value = true;
         } else {
             return_value = false;
