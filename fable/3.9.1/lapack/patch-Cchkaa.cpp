--- lin/common/Cchkaa.cpp_	2026-01-25 13:32:01.520723708 +0900
+++ lin/common/Cchkaa.cpp	2026-01-25 13:32:01.526723872 +0900
@@ -43,8 +43,10 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_zchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+#include <memory>
+
+void Cchkaa(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
@@ -54,10 +56,13 @@
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
@@ -93,12 +98,18 @@
     INTEGER ntypes = 0;
     bool dotype[matmax];
     const INTEGER kdmax = nmax + (nmax + 1) / 4;
-    COMPLEX a[(kdmax + 1) * nmax * 7];
-    COMPLEX b[nmax * maxrhs * 4];
-    COMPLEX work[nmax * nmax + maxrhs + 10];
-    REAL rwork[150 + 2 * maxrhs];
-    INTEGER iwork[25 * nmax];
-    REAL s[2 * nmax];
+    auto a_storage = std::make_unique<COMPLEX[]>((kdmax + 1) * nmax * 7);
+    auto b_storage = std::make_unique<COMPLEX[]>(nmax * maxrhs * 4);
+    auto work_storage = std::make_unique<COMPLEX[]>(nmax * nmax + maxrhs + 10);
+    auto rwork_storage = std::make_unique<REAL[]>(150 + 2 * maxrhs);
+    auto iwork_storage = std::make_unique<INTEGER[]>(25 * nmax);
+    auto s_storage = std::make_unique<REAL[]>(2 * nmax);
+    COMPLEX *a = a_storage.get();
+    COMPLEX *b = b_storage.get();
+    COMPLEX *work = work_storage.get();
+    REAL *rwork = rwork_storage.get();
+    INTEGER *iwork = iwork_storage.get();
+    REAL *s = s_storage.get();
     INTEGER la = 0;
     INTEGER lafac = 0;
     INTEGER piv[nmax];
@@ -112,8 +123,9 @@
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
@@ -133,8 +145,8 @@
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
@@ -1082,4 +1094,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkaa); }
+int main(int argc, char const *argv[]) { Cchkaa(); }
