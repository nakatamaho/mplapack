--- Cchkaa.cpp_	2026-01-30 06:56:30.290306525 +0900
+++ Cchkaa.cpp	2026-01-30 06:56:44.249998327 +0900
@@ -44,8 +44,8 @@
 #include <mplapack_lin.h>
 #include <memory>
 
-void program_zchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkaa(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
@@ -55,10 +55,13 @@
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
@@ -119,8 +123,9 @@
     static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
     static const char *format_9996 = "(' Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
     static const char *format_9995 = "(' Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
-    static const char *format_9994 = "(' Tests of the COMPLEX*16 LAPACK routines ',/,' LAPACK VERSION ',i1,'.',"
-                                     "i1,'.',i1,/,/,' The following parameter values will be used:')";
+    static const char *format_9994 = "(' Tests of the complex multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
+                                     "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "' The following parameter values will be used:')";
     static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
     static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                      "f8.2,/)";
@@ -140,8 +145,8 @@
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
@@ -440,7 +445,7 @@
     //
     // Check first character for correct precision.
     //
-    if (!Mlsame(c1.elems, "Zomplex precision")) {
+    if (!Mlsame(c1.elems, "Zomplex precision") && !Mlsame(c1.elems, "C")) {
         write(nout, format_9990), path;
         //
     } else if (nmats <= 0) {
@@ -1085,7 +1085,7 @@ statement_140:
     cmn.io.close(nin);
     s2 = dsecnd();
     write(nout, format_9998);
-    write(nout, format_9997), s2 - s1;
+    write(nout, format_9997), cast2double(s2 - s1);
     //
     // End of Cchkaa
     //
@@ -1089,4 +1094,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkaa); }
+int main(int argc, char const *argv[]) { Cchkaa(); }
diff --git a/mplapack/test/lin/common/Cchkaa.cpp b/mplapack/test/lin/common/Cchkaa.cpp
