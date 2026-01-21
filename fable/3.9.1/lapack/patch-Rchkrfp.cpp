--- Rchkrfp.cpp_	2026-01-21 18:29:14.733032980 +0900
+++ Rchkrfp.cpp	2026-01-21 18:29:14.740033112 +0900
@@ -43,7 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_dchkrfp(int argc, char const *argv[]) {
+void program_Rchkrfp(int argc, char const *argv[]) {
     common cmn(argc, argv);
     common_read read(cmn);
     common_write write(cmn);
@@ -60,18 +60,21 @@
     // Read a dummy line.
     //
     const INTEGER nin = 5;
+    const INTEGER nout = 6;
     read(nin, star);
     //
     // Report LAPACK version tag (e.g. LAPACK-3.2.0)
     //
-    INTEGER vers_major = 0;
-    INTEGER vers_minor = 0;
-    INTEGER vers_patch = 0;
-    ilaver(vers_major, vers_minor, vers_patch);
-    const INTEGER nout = 6;
-    write(nout, "(/,' Tests of the DOUBLE PRECISION LAPACK RFP routines ',/,"
-                "' LAPACK VERSION ',i1,'.',i1,'.',i1,/,/,"
-                "' The following parameter values will be used:')"),
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, "(' Tests of the Multiple precision version of LAPACK RFP routines',i1,'.',i1,'.',i1,/, "
+	  "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
+        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of N
@@ -279,4 +282,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkrfp); }
+int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_Rchkrfp); }
