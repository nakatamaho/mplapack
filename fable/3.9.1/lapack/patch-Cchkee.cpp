--- eig/common/Cchkee.cpp_	2026-01-25 13:49:45.074727946 +0900
+++ eig/common/Cchkee.cpp	2026-01-25 13:49:45.081728106 +0900
@@ -43,19 +43,17 @@
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
-void program_zchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+#include <memory>
+
+void Cchkee(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static INTEGER ioldsd[4] = {0, 0, 0, 1};
     static fem::str<10> intstr = "0123456789";
     const INTEGER nmax = 132;
     const INTEGER need = 14;
-    COMPLEX a[nmax * nmax * need];
-    COMPLEX b[nmax * nmax * 5];
     const INTEGER ncmax = 20;
-    COMPLEX c[ncmax * ncmax * ncmax * ncmax];
-    COMPLEX dc[nmax * 6];
     REAL s1 = 0.0;
     bool fatal = false;
     const INTEGER nout = 6;
@@ -87,9 +85,12 @@
     bool zgk = false;
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
@@ -128,23 +129,32 @@
     INTEGER maxtyp = 0;
     bool dotype[maxt];
     const INTEGER lwork = nmax * (5 * nmax + 20);
-    COMPLEX work[lwork];
-    REAL rwork[lwork];
     const INTEGER liwork = nmax * (nmax + 20);
     INTEGER iwork[liwork];
     bool logwrk[nmax];
-    REAL result[500];
     INTEGER info = 0;
     REAL dr[nmax * 12];
     INTEGER nrhs = 0;
     bool tstdif = false;
     REAL thrshn = 0.0;
-    REAL s[nmax * nmax];
-    COMPLEX x[5 * nmax];
-    COMPLEX taua[nmax];
-    COMPLEX taub[nmax];
-    REAL alpha[nmax];
-    REAL beta[nmax];
+    auto work_storage = std::make_unique<COMPLEX[]>(lwork);
+    auto rwork_storage = std::make_unique<REAL[]>(lwork);
+    auto result_storage = std::make_unique<REAL[]>(500);
+    auto s_storage = std::make_unique<REAL[]>(nmax * nmax);
+    auto x_storage = std::make_unique<COMPLEX[]>(5 * nmax);
+    auto taua_storage = std::make_unique<COMPLEX[]>(nmax);
+    auto taub_storage = std::make_unique<COMPLEX[]>(nmax);
+    auto alpha_storage = std::make_unique<REAL[]>(nmax);
+    auto beta_storage = std::make_unique<REAL[]>(nmax);
+    COMPLEX *work = work_storage.get();
+    REAL *rwork = rwork_storage.get();
+    REAL *result = result_storage.get();
+    REAL *s = s_storage.get();
+    COMPLEX *x = x_storage.get();
+    COMPLEX *taua = taua_storage.get();
+    COMPLEX *taub = taub_storage.get();
+    REAL *alpha = alpha_storage.get();
+    REAL *beta = beta_storage.get();
     REAL s2 = 0.0;
     INTEGER lda = nmax * nmax;
     INTEGER ldb = nmax * nmax;
@@ -184,7 +194,9 @@
     static const char *format_9974 = "(' Tests of Chbtrd',/,' (reduction of a Hermitian band ',"
                                      "'matrix to real tridiagonal form)')";
     static const char *format_9973 = "(/,1x,71('-'))";
-    static const char *format_9972 = "(/,' LAPACK VERSION ',i1,'.',i1,'.',i1)";
+    static const char *format_9972 = "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
+                                     "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "'The following parameter values will be used:')";
     static const char *format_9971 = "(/,' Tests of the Generalized Linear Regression Model ','routines')";
     static const char *format_9970 = "(/,' Tests of the Generalized QR and RQ routines')";
     static const char *format_9969 = "(/,' Tests of the Generalized Singular Value',' Decomposition routines')";
@@ -205,10 +217,19 @@
     static const char *format_9960 = "(/,' Tests of the CS Decomposition routines')";
     //
     //
-    a = 0.0;
-    b = 0.0;
-    c = 0.0;
-    dc = 0.0;
+    auto a_storage = std::make_unique<COMPLEX[]>(nmax * nmax * need);
+    auto b_storage = std::make_unique<COMPLEX[]>(nmax * nmax * 5);
+    auto c_storage = std::make_unique<COMPLEX[]>(ncmax * ncmax * ncmax * ncmax);
+    auto dc_storage = std::make_unique<COMPLEX[]>(nmax * 6);
+    COMPLEX *a = a_storage.get();
+    COMPLEX *b = b_storage.get();
+    COMPLEX *c = c_storage.get();
+    COMPLEX *dc = dc_storage.get();
+    const COMPLEX zero = COMPLEX(0.0);
+    std::fill_n(a, nmax * nmax * need, zero);
+    std::fill_n(b, nmax * nmax * 5, zero);
+    std::fill_n(c, ncmax * ncmax * ncmax * ncmax, zero);
+    std::fill_n(dc, nmax * 6, zero);
     s1 = dsecnd();
     fatal = false;
     nunit = nout;
@@ -329,8 +350,8 @@
         write(nout, format_9992), path;
         goto statement_380;
     }
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, format_9972), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9972), mplapack_vers_major, mplapack_vers_minor, mplapck_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     write(nout, format_9984);
     //
     // Read the number of values of M, P, and N.
@@ -1054,7 +1075,7 @@
                     iseed[k - 1] = ioldsd[k - 1];
                 }
             }
-            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], fem::max(11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
             Cchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &dc[(3 - 1) * nmax], work, lwork, rwork, iwork, logwrk, result, info);
             if (info != 0) {
                 write(nout, format_9980), "Cchkhs", info;
@@ -1590,4 +1611,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkee); }
+int main(int argc, char const *argv[]) { Cchkee(); }
