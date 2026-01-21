--- Rchkaa.cpp_	2026-01-21 18:28:37.385333136 +0900
+++ Rchkaa.cpp	2026-01-21 18:28:37.389333210 +0900
@@ -43,7 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_dchkaa(int argc, char const *argv[]) {
+void program_Rchkaa(int argc, char const *argv[]) {
     common cmn(argc, argv);
     common_read read(cmn);
     common_write write(cmn);
@@ -54,10 +54,13 @@
     INTEGER lda = 0;
     bool fatal = false;
     const INTEGER nin = 5;
-    INTEGER vers_major = 0;
-    INTEGER vers_minor = 0;
-    INTEGER vers_patch = 0;
     const INTEGER nout = 6;
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
     INTEGER nm = 0;
     const INTEGER maxin = 12;
     INTEGER mval[maxin];
@@ -125,10 +128,10 @@
     //
     // Report values of parameters.
     //
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, "(' Tests of the DOUBLE PRECISION LAPACK routines ',/,' LAPACK VERSION ',"
-                "i1,'.',i1,'.',i1,/,/,' The following parameter values will be used:')"),
-        vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
+                "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
+        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of M
     //
@@ -931,4 +934,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkaa); }
+int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_Rchkaa); }
