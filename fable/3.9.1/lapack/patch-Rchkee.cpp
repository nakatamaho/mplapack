--- Rchkee.cpp
+++ Rchkee.cpp
@@ -43,8 +43,8 @@ using fem::common;
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
-void program_dchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Rchkee() {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static INTEGER ioldsd[4] = {0, 0, 0, 1};
@@ -87,9 +87,12 @@ void program_dchkee(int argc, char const *argv[]) {
     bool dgk = false;
     REAL thresh = 0.0;
     bool tsterr = false;
-    INTEGER vers_major = 0;
-    INTEGER vers_minor = 0;
-    INTEGER vers_patch = 0;
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
     INTEGER nn = 0;
     const INTEGER maxin = 20;
     INTEGER mval[maxin];
@@ -296,8 +299,10 @@ statement_10:
         write(nout, format_9992), path;
         goto statement_10;
     }
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, "(/,' LAPACK VERSION ',i1,'.',i1,'.',i1)"), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
+                "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
+        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     write(nout, "(/,' The following parameter values will be used:')");
     //
     // Read the number of values of M, P, and N.
@@ -1563,4 +1568,4 @@ statement_380:
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkee); }
+int main(int argc, char const *argv[]) { Rchkee(); }
