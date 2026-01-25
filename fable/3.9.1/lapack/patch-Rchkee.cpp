--- eig/common/Rchkee.cpp_	2026-01-25 13:32:01.498723107 +0900
+++ eig/common/Rchkee.cpp	2026-01-25 13:32:01.506723325 +0900
@@ -42,20 +42,18 @@
 
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
+#include <mplapack_debug.h>
 
-void program_dchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+#include <memory>
+void Rchkee(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static INTEGER ioldsd[4] = {0, 0, 0, 1};
     static fem::str<10> intstr = "0123456789";
     const INTEGER nmax = 132;
-    const INTEGER need = 14;
-    REAL a[nmax * nmax * need];
-    REAL b[nmax * nmax * 5];
     const INTEGER ncmax = 20;
-    REAL c[ncmax * ncmax * ncmax * ncmax];
-    REAL d[nmax * 12];
+    const INTEGER need = 14;
     REAL s1 = 0.0;
     bool fatal = false;
     const INTEGER nout = 6;
@@ -87,9 +85,12 @@
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
@@ -128,18 +129,23 @@
     INTEGER maxtyp = 0;
     bool dotype[maxt];
     const INTEGER lwork = nmax * (5 * nmax + 5) + 1;
-    REAL work[lwork];
     const INTEGER liwork = nmax * (5 * nmax + 20);
     INTEGER iwork[liwork];
     bool logwrk[nmax];
-    REAL result[500];
     INTEGER info = 0;
     INTEGER nrhs = 0;
     bool tstdif = false;
     REAL thrshn = 0.0;
-    REAL x[5 * nmax];
-    REAL taua[nmax];
-    REAL taub[nmax];
+    auto result_storage = std::make_unique<REAL[]>(500);
+    auto work_storage = std::make_unique<REAL[]>(lwork);
+    auto x_storage = std::make_unique<REAL[]>(5 * nmax);
+    auto taua_storage = std::make_unique<REAL[]>(nmax);
+    auto taub_storage = std::make_unique<REAL[]>(nmax);
+    REAL *result = result_storage.get();
+    REAL *work = work_storage.get();
+    REAL *x = x_storage.get();
+    REAL *taua = taua_storage.get();
+    REAL *taub = taub_storage.get();
     REAL s2 = 0.0;
     INTEGER lda = nmax * nmax;
     INTEGER ldc = ncmax * ncmax;
@@ -180,7 +186,9 @@
     static const char *format_9974 = "(' Tests of Rsbtrd',/,' (reduction of a symmetric band ',"
                                      "'matrix to tridiagonal form)')";
     static const char *format_9973 = "(/,1x,71('-'))";
-    static const char *format_9972 = "(/,' LAPACK VERSION ',i1,'.',i1,'.',i1)";
+    static const char *format_9972 = "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
+                                     "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "'The following parameter values will be used:')";
     static const char *format_9971 = "(/,' Tests of the Generalized Linear Regression Model ','routines')";
     static const char *format_9970 = "(/,' Tests of the Generalized QR and RQ routines')";
     static const char *format_9969 = "(/,' Tests of the Generalized Singular Value',' Decomposition routines')";
@@ -201,13 +209,16 @@
     static const char *format_9960 = "(/,' Tests of the CS Decomposition routines')";
     //
     //
-    a = 0.0;
-    b = 0.0;
-    c = 0.0;
-    d = 0.0;
+    auto a_storage = std::make_unique<REAL[]>(nmax * nmax * need);
+    auto b_storage = std::make_unique<REAL[]>(nmax * nmax * 5);
+    auto c_storage = std::make_unique<REAL[]>(ncmax * ncmax * ncmax * ncmax);
+    auto d_storage = std::make_unique<REAL[]>(nmax * 12);
+    REAL *a = a_storage.get();
+    REAL *b = b_storage.get();
+    REAL *c = c_storage.get();
+    REAL *d = d_storage.get();
     s1 = dsecnd();
     fatal = false;
-    nunit = nout;
 //
 // Return to here to read multiple sets of data
 //
@@ -329,8 +340,8 @@
         write(nout, format_9992), path;
         goto statement_10;
     }
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, format_9972), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9994), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     write(nout, format_9984);
     //
     // Read the number of values of M, P, and N.
@@ -1588,4 +1599,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkee); }
+int main(int argc, char const *argv[]) { Rchkee(); }
