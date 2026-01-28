--- Rchkaa.cpp_	2026-01-28 09:31:30.425494060 +0900
+++ Rchkaa.cpp	2026-01-28 09:31:30.431494190 +0900
@@ -43,8 +43,10 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_dchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+#include <memory>
+
+void Rchkaa(void) {
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
@@ -93,16 +98,22 @@
     INTEGER ntypes = 0;
     bool dotype[matmax];
     const INTEGER kdmax = nmax + (nmax + 1) / 4;
-    REAL a[(kdmax + 1) * nmax * 7];
-    REAL b[nmax * maxrhs * 4];
-    REAL work[nmax * 3 * nmax + maxrhs + 30];
-    REAL rwork[5 * nmax + 2 * maxrhs];
+    auto a_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (kdmax + 1) * nmax * 7));
+    REAL *a = a_storage.get();
+    auto b_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs * 4));
+    REAL *b = b_storage.get();
+    auto work_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * 3 * nmax + maxrhs + 30));
+    REAL *work = work_storage.get();
+    auto rwork_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 5 * nmax + 2 * maxrhs));
+    REAL *rwork = rwork_storage.get();
+    auto s_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 2 * nmax));
+    REAL *s = s_storage.get();
+    auto e_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *e = e_storage.get();
     INTEGER iwork[25 * nmax];
-    REAL s[2 * nmax];
     INTEGER la = 0;
     INTEGER lafac = 0;
     INTEGER piv[nmax];
-    REAL e[nmax];
     REAL s2 = 0.0;
     INTEGER ldaw = (kdmax + 1) * nmax;
     INTEGER ldb = nmax * maxrhs;
@@ -112,8 +123,9 @@
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
@@ -434,7 +446,7 @@
     //
     // Check first character for correct precision.
     //
-    if (!Mlsame(c1.elems, "Double precision")) {
+    if (!Mlsame(c1.elems, "Double precision") && !Mlsame(c1.elems, "R")) {
         write(nout, format_9990), path;
         //
     } else if (nmats <= 0) {
@@ -935,4 +947,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkaa); }
+int main(int argc, char const *argv[]) { Rchkaa(); }
