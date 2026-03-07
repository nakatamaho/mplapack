--- Rchkaa.cpp_	2026-01-30 07:01:53.559480665 +0900
+++ Rchkaa.cpp	2026-01-30 07:01:59.392463469 +0900
@@ -42,10 +42,11 @@
 
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
+
 #include <memory>
 
-void program_dchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Rchkaa(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
@@ -56,10 +56,13 @@
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
@@ -120,8 +123,9 @@
     static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
     static const char *format_9996 = "(' Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
     static const char *format_9995 = "(' Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
-    static const char *format_9994 = "(' Tests of the DOUBLE PRECISION LAPACK routines ',/,' LAPACK VERSION ',"
-                                     "i1,'.',i1,'.',i1,/,/,' The following parameter values will be used:')";
+    static const char *format_9994 = "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
+                                     "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "' The following parameter values will be used:')";
     static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
     static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                      "f8.2,/)";
@@ -141,8 +145,8 @@
     //
     // Report values of parameters.
     //
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, format_9994), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9994), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of M
     //
@@ -442,7 +446,7 @@
     //
     // Check first character for correct precision.
     //
-    if (!Mlsame(c1.elems, "Double precision")) {
+    if (!Mlsame(c1.elems, "Double precision") && !Mlsame(c1.elems, "R")) {
         write(nout, format_9990), path;
         //
     } else if (nmats <= 0) {
@@ -939,7 +939,7 @@ statement_140:
     cmn.io.close(nin);
     s2 = dsecnd();
     write(nout, format_9998);
-    write(nout, format_9997), s2 - s1;
+    write(nout, format_9997), cast2double(s2 - s1);
     //
     // End of Rchkaa
     //
@@ -943,4 +947,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkaa); }
+int main(int argc, char const *argv[]) { Rchkaa(); }
