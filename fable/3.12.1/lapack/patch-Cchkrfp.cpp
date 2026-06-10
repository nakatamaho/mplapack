--- Cchkrfp.cpp_	2026-01-28 09:31:30.411493758 +0900
+++ Cchkrfp.cpp	2026-01-28 09:31:30.417493887 +0900
@@ -44,8 +44,8 @@
 #include <mplapack_lin.h>
 #include <memory>
 
-void program_zchkrfp(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkrfp(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     //
@@ -54,8 +54,9 @@
     static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
     static const char *format_9996 = "(' !! Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
     static const char *format_9995 = "(' !! Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
-    static const char *format_9994 = "(/,' Tests of the COMPLEX*16 LAPACK RFP routines ',/,' LAPACK VERSION ',"
-                                     "i1,'.',i1,'.',i1,/,/,' The following parameter values will be used:')";
+    static const char *format_9994 = "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
+                                     "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "'The following parameter values will be used:')";
     static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
     static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                      "f8.2,/)";
@@ -68,16 +69,19 @@
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
-    write(nout, format_9994), vers_major, vers_minor, vers_patch;
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9994), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of N
     //
@@ -294,7 +294,7 @@ void Cchkrfp(void) {
     cmn.io.close(nin);
     REAL s2 = dsecnd();
     write(nout, format_9998);
-    write(nout, format_9997), s2 - s1;
+    write(nout, format_9997), cast2double(s2 - s1);
     //
     // End of Cchkrfp
     //
@@ -302,4 +306,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkrfp); }
+int main(int argc, char const *argv[]) { Cchkrfp(); }
