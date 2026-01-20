--- Cchkaa.cpp_	2026-01-20 13:35:36.236672493 +0900
+++ Cchkaa.cpp	2026-01-20 13:35:41.691787391 +0900
@@ -43,8 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_zchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkaa(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
@@ -54,9 +54,12 @@
     INTEGER lda = 0;
     bool fatal = false;
     const INTEGER nin = 5;
-    INTEGER vers_major = 0;
-    INTEGER vers_minor = 0;
-    INTEGER vers_patch = 0;
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
     const INTEGER nout = 6;
     INTEGER nm = 0;
     const INTEGER maxin = 12;
@@ -125,10 +125,10 @@
     //
     // Report values of parameters.
     //
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, "(' Tests of the COMPLEX*16 LAPACK routines ',/,' LAPACK VERSION ',i1,'.',"
-                "i1,'.',i1,/,/,' The following parameter values will be used:')"),
-        vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
+                "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
+        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of M
     //
@@ -1078,4 +1078,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkaa); }
+int main(int argc, char const *argv[]) { Cchkaa(); }
