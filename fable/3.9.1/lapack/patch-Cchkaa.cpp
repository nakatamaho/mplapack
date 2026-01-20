--- Cchkaa.cpp_	2026-01-20 13:06:35.630160802 +0900
+++ Cchkaa.cpp	2026-01-20 13:06:40.900284867 +0900
@@ -43,8 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_zchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkaa(void) {
     common_read read(cmn);
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
@@ -125,10 +124,10 @@
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
@@ -1078,4 +1077,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkaa); }
+int main(int argc, char const *argv[]) { return Cchkaa(); }
