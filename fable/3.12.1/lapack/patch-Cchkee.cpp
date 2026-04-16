--- Cchkee.cpp	2026-02-17 19:46:40.127204043 +0900
+++ Cchkee.cpp	2026-02-17 19:40:45.260876385 +0900
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
@@ -99,10 +98,14 @@
     bool zgl = false;
     bool zgk = false;
     REAL thresh = 0.0;
+    REAL thresh_org = 0.0;
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
@@ -194,7 +197,9 @@
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
@@ -214,47 +219,24 @@
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
@@ -375,8 +357,8 @@
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
@@ -1100,7 +1082,17 @@
                     iseed[k - 1] = ioldsd[k - 1];
                 }
             }
-            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], fem::max(11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 5.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Cchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &dc[(3 - 1) * nmax], work, lwork, rwork, iwork, logwrk, result, info);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
             if (info != 0) {
                 write(nout, format_9980), "Cchkhs", info;
            }
@@ -1185,4 +1221,14 @@
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                thresh_org = thresh;
+                thresh = thresh * 2.0;
+                printf("Warning! Threshold has been lifted to: ");
+                printnum_short(thresh);
+                printf(" for GMP ZSG\n");
+#endif
                 Cdrvsg2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], nmax, &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &a[(7 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+                thresh = thresh_org;
+#endif
                 if (info != 0) {
                     write(nout, format_9980), "Cdrvsg", info;
                 }
@@ -1234,6 +1226,13 @@
                 }
             }
             write(nout, format_9995), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], nrhs;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 4.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             if (tstchk) {
                 Cchkbd(nn, mval, nval, maxtyp, dotype, nrhs, iseed, thresh, &a[0], nmax, &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], nmax, &a[(7 - 1) * lda], &a[(8 - 1) * lda], work, lwork, rwork, nout, info);
                 if (info != 0) {
@@ -1243,6 +1242,9 @@
             if (tstdrv) {
                 Cdrvbd(nn, mval, nval, maxtyp, dotype, iseed, thresh, &a[0], nmax, &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], work, lwork, rwork, iwork, nout, info);
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         //
     } else if (Mlsamen(3, c3.elems, "ZEV")) {
@@ -1260,11 +1262,21 @@
             if (tsterr) {
                 Cerred(c3, nout);
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 5.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Cdrvev(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, result, work, lwork, rwork, iwork, info);
             if (info != 0) {
                 write(nout, format_9980), "Cgeev", info;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1284,11 +1296,21 @@
             if (tsterr) {
                 Cerred(c3, nout);
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 5.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Cdrves(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(4 - 1) * lda], nmax, result, work, lwork, rwork, iwork, logwrk, info);
             if (info != 0) {
                 write(nout, format_9980), "Cgees", info;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1308,11 +1330,21 @@
             if (tsterr) {
                 Cerred(c3, nout);
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 4.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Cdrvvx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(6 - 1) * nmax], &dr[(7 - 1) * nmax], &dr[(8 - 1) * nmax], result, work, lwork, rwork, info);
             if (info != 0) {
                 write(nout, format_9980), "Cgeevx", info;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1332,11 +1364,21 @@
             if (tsterr) {
                 Cerred(c3, nout);
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 5.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Cdrvsx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &dc[(3 - 1) * nmax], &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], result, work, lwork, rwork, logwrk, info);
             if (info != 0) {
                 write(nout, format_9980), "Cgeesx", info;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1401,6 +1443,12 @@
                 Cerrgg(c3, nout);
             }
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh * 5.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Cdrges(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(7 - 1) * lda], nmax, &a[(8 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], work, lwork, rwork, result, logwrk, info);
             //
             if (info != 0) {
@@ -1414,6 +1462,9 @@
             if (info != 0) {
                 write(nout, format_9980), "Cdrges3", info;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1459,6 +1510,13 @@
                 Cerrgg(c3, nout);
             }
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh_org = thresh;
+            thresh = thresh * 4.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Cdrgev(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(7 - 1) * lda], nmax, &a[(8 - 1) * lda], &a[(9 - 1) * lda], nmax, &dc[0], &dc[(2 - 1) * nmax], &dc[(3 - 1) * nmax], &dc[(4 - 1) * nmax], work, lwork, rwork, result, info);
             if (info != 0) {
                 write(nout, format_9980), "Cdrgev", info;
@@ -1471,6 +1528,9 @@
             if (info != 0) {
                 write(nout, format_9980), "Cdrgev3", info;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1630,10 +1690,10 @@
 statement_380:
     write(nout, format_9994);
     s2 = dsecnd();
-    write(nout, format_9993), s2 - s1;
+    write(nout, format_9993), cast2double(s2 - s1);
     //
     // End of Cchkee
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkee); }
+int main(int argc, char const *argv[]) { Cchkee(); }
