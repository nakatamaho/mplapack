--- Cchkrfp.cpp_	2026-01-21 19:57:58.924315069 +0900
+++ Cchkrfp.cpp	2026-01-21 19:57:58.931315193 +0900
@@ -43,8 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_zchkrfp(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkrfp(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static const char *format_9991 = "(' Relative machine ',a,' is taken to be',d16.6)";
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
-    write(nout, "(/,' Tests of the COMPLEX*16 LAPACK RFP routines ',/,' LAPACK VERSION ',"
-                "i1,'.',i1,'.',i1,/,/,' The following parameter values will be used:')"),
-        vers_major, vers_minor, vers_patch;
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
+	  "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
+        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of N
     //
@@ -279,4 +282,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkrfp); }
+int main(int argc, char const *argv[]) { Cchkrfp(); }
