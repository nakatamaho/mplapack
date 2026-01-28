--- Cchkaa.cpp_	2026-01-28 09:31:30.397493455 +0900
+++ Cchkaa.cpp	2026-01-28 09:31:30.403493585 +0900
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
@@ -93,16 +98,23 @@
     INTEGER ntypes = 0;
     bool dotype[matmax];
     const INTEGER kdmax = nmax + (nmax + 1) / 4;
-    COMPLEX a[(kdmax + 1) * nmax * 7];
-    COMPLEX b[nmax * maxrhs * 4];
-    COMPLEX work[nmax * nmax + maxrhs + 10];
-    REAL rwork[150 + 2 * maxrhs];
-    INTEGER iwork[25 * nmax];
-    REAL s[2 * nmax];
+    auto a_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, (kdmax + 1) * nmax * 7));
+    COMPLEX *a = a_storage.get();
+    auto b_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * maxrhs * 4));
+    COMPLEX *b = b_storage.get();
+    auto work_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * nmax + maxrhs + 10));
+    COMPLEX *work = work_storage.get();
+    auto rwork_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 150 + 2 * maxrhs));
+    REAL *rwork = rwork_storage.get();
+    auto iwork_storage = std::make_unique<INTEGER[]>(std::max<INTEGER>(1, 25 * nmax));
+    INTEGER *iwork = iwork_storage.get();
+    auto s_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 2 * nmax));
+    REAL *s = s_storage.get();
+    auto e_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax));
+    COMPLEX *e = e_storage.get();
     INTEGER la = 0;
     INTEGER lafac = 0;
     INTEGER piv[nmax];
-    COMPLEX e[nmax];
     REAL s2 = 0.0;
     INTEGER ldaw = (kdmax + 1) * nmax;
     INTEGER ldb = nmax * maxrhs;
@@ -112,8 +124,9 @@
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
@@ -133,8 +146,8 @@
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
@@ -433,7 +446,7 @@
     //
     // Check first character for correct precision.
     //
-    if (!Mlsame(c1.elems, "Zomplex precision")) {
+    if (!Mlsame(c1.elems, "Zomplex precision") && !Mlsame(c1.elems, "C")) {
         write(nout, format_9990), path;
         //
     } else if (nmats <= 0) {
@@ -1082,4 +1095,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkaa); }
+int main(int argc, char const *argv[]) { Cchkaa(); }
