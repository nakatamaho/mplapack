--- Cchkee.cpp_	2026-01-28 11:33:05.861116874 +0900
+++ Cchkee.cpp	2026-01-28 11:33:05.865116941 +0900
@@ -44,13 +44,12 @@
 #include <mplapack_eig.h>
 #include <memory>
 
-void program_zchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkee(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static INTEGER ioldsd[4] = {0, 0, 0, 1};
     static fem::str<10> intstr = "0123456789";
-    INTEGER allocatestatus = 0;
     const INTEGER nmax = 132;
     std::unique_ptr<REAL[]> s_storage;
     REAL *s = nullptr;
@@ -100,9 +99,12 @@
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
@@ -200,7 +202,9 @@
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
@@ -220,47 +224,24 @@
                                      "', INWIN =',i4,', INIBL =',i4,', ISHFTS =',i4,', IACC22 =',i4)";
     static const char *format_9960 = "(/,' Tests of the CS Decomposition routines')";
     //
-    allocatestatus = 0;
     s_storage = std::make_unique<REAL[]>(max((INTEGER)1, (nmax * nmax)));
     s = s_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     a_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (nmax * nmax) * need));
     a = a_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     b_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (nmax * nmax) * 5));
     b = b_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     c_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (ncmax * ncmax) * (ncmax * ncmax)));
     c = c_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
     rwork = rwork_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lwork));
     work = work_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
     //
-    a = 0.0;
-    b = 0.0;
-    c = 0.0;
-    dc = 0.0;
+    const COMPLEX zero = COMPLEX(0.0);
+    std::fill_n(a, nmax * nmax * need, zero);
+    std::fill_n(b, nmax * nmax * 5, zero);
+    std::fill_n(c, ncmax * ncmax * ncmax * ncmax, zero);
+    std::fill_n(dc, nmax * 6, zero);
     s1 = dsecnd();
     fatal = false;
     nunit = nout;
@@ -381,8 +362,8 @@
         write(nout, format_9992), path;
         goto statement_380;
     }
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, format_9972), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9972), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     write(nout, format_9984);
     //
     // Read the number of values of M, P, and N.
@@ -1106,7 +1087,7 @@
                     iseed[k - 1] = ioldsd[k - 1];
                 }
             }
-            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], fem::max(11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
             Cchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &dc[(3 - 1) * nmax], work, lwork, rwork, iwork, logwrk, result, info);
             if (info != 0) {
                 write(nout, format_9980), "Cchkhs", info;
@@ -1611,7 +1611,7 @@ statement_190:
 statement_380:
     write(nout, format_9994);
     s2 = dsecnd();
-    write(nout, format_9993), s2 - s1;
+    write(nout, format_9993), cast2double(s2 - s1);
     //
     // End of Cchkee
     //
@@ -1642,4 +1623,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkee); }
+int main(int argc, char const *argv[]) { Cchkee(); }
