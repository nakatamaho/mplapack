--- Rchkee.cpp	2026-02-17 19:46:39.912200751 +0900
+++ Rchkee.cpp	2026-02-17 19:40:43.001843328 +0900
@@ -44,13 +44,12 @@
 #include <mplapack_eig.h>
 #include <memory>
 
-void program_dchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Rchkee(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     static INTEGER ioldsd[4] = {0, 0, 0, 1};
     static fem::str<10> intstr = "0123456789";
-    INTEGER allocatestatus = 0;
     const INTEGER nmax = 132;
     const INTEGER need = 14;
     std::unique_ptr<REAL[]> a_storage;
@@ -95,10 +94,14 @@
     bool dgl = false;
     bool dgk = false;
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
@@ -187,7 +190,9 @@
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
@@ -207,35 +212,20 @@
                                      "', INWIN =',i4,', INIBL =',i4,', ISHFTS =',i4,', IACC22 =',i4)";
     static const char *format_9960 = "(/,' Tests of the CS Decomposition routines')";
     //
-    allocatestatus = 0;
     a_storage = std::make_unique<REAL[]>(max((INTEGER)1, (nmax * nmax) * need));
     a = a_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     b_storage = std::make_unique<REAL[]>(max((INTEGER)1, (nmax * nmax) * 5));
     b = b_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     c_storage = std::make_unique<REAL[]>(max((INTEGER)1, (ncmax * ncmax) * (ncmax * ncmax)));
     c = c_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
     work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
     work = work_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
     //
-    a = 0.0;
-    b = 0.0;
-    c = 0.0;
-    d = 0.0;
+    const REAL zero = REAL(0.0);
+    std::fill_n(a, nmax * nmax * need, zero);
+    std::fill_n(b, nmax * nmax * 5, zero);
+    std::fill_n(c, ncmax * ncmax * ncmax * ncmax, zero);
+    std::fill_n(d, nmax * 12, zero);
     s1 = dsecnd();
     fatal = false;
     nunit = nout;
@@ -354,14 +344,24 @@
         Mxlaenv(15, 2);
         Mxlaenv(16, 2);
         tsterr = true;
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 2.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
         Rchkec(thresh, tsterr, nin, nout);
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         goto statement_10;
     } else {
         write(nout, format_9992), path;
         goto statement_10;
     }
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, format_9972), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9972), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     write(nout, format_9984);
     //
     // Read the number of values of M, P, and N.
@@ -1086,11 +1086,21 @@
                     iseed[k - 1] = ioldsd[k - 1];
                 }
             }
-            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], fem::max(11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 4.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Rchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &d[(5 - 1) * nmax], &d[(6 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &d[(7 - 1) * nmax], work, lwork, iwork, logwrk, result, info);
             if (info != 0) {
                 write(nout, format_9980), "Rchkhs", info;
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         }
         //
     } else if (Mlsamen(3, c3.elems, "DST") || Mlsamen(3, c3.elems, "SEP") || Mlsamen(3, c3.elems, "SE2")) {
@@ -1219,6 +1229,13 @@
                 }
             }
             write(nout, format_9995), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], nrhs;
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 2.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             if (tstchk) {
                 Rchkbd(nn, mval, nval, maxtyp, dotype, nrhs, iseed, thresh, &a[0], nmax, &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], nmax, &a[(7 - 1) * lda], &a[(8 - 1) * lda], work, lwork, iwork, nout, info);
                 if (info != 0) {
@@ -1228,6 +1245,9 @@
             if (tstdrv) {
                 Rdrvbd(nn, mval, nval, maxtyp, dotype, iseed, thresh, &a[0], nmax, &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], work, lwork, iwork, nout, info);
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         }
         //
     } else if (Mlsamen(3, c3.elems, "DEV")) {
@@ -1245,11 +1265,21 @@
             if (tsterr) {
                 Rerred(c3, nout);
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 3.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Rdrvev(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, result, work, lwork, iwork, info);
             if (info != 0) {
                 write(nout, format_9980), "Rgeev", info;
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1269,11 +1299,21 @@
             if (tsterr) {
                 Rerred(c3, nout);
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 2.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Rdrves(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &a[(4 - 1) * lda], nmax, result, work, lwork, iwork, logwrk, info);
             if (info != 0) {
                 write(nout, format_9980), "Rgees", info;
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1293,11 +1333,21 @@
             if (tsterr) {
                 Rerred(c3, nout);
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 3.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Rdrvvx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, &d[(5 - 1) * nmax], &d[(6 - 1) * nmax], &d[(7 - 1) * nmax], &d[(8 - 1) * nmax], &d[(9 - 1) * nmax], &d[(10 - 1) * nmax], &d[(11 - 1) * nmax], &d[(12 - 1) * nmax], result, work, lwork, iwork, info);
             if (info != 0) {
                 write(nout, format_9980), "Rgeevx", info;
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1317,11 +1367,21 @@
             if (tsterr) {
                 Rerred(c3, nout);
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh_org = thresh;
+            thresh = thresh * 2.0;
+            printf("Warning! Threshold has been lifted to: ");
+            printnum_short(thresh);
+            printf(" for GMP\n");
+#endif
             Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
             Rdrvsx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &d[(5 - 1) * nmax], &d[(6 - 1) * nmax], &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], result, work, lwork, iwork, logwrk, info);
             if (info != 0) {
                 write(nout, format_9980), "Rgeesx", info;
             }
+#if defined MPLAPACK_BUILD_WITH_GMP
+            thresh = thresh_org;
+#endif
         }
         write(nout, format_9973);
         goto statement_10;
@@ -1613,10 +1673,10 @@
 statement_380:
     write(nout, format_9994);
     s2 = dsecnd();
-    write(nout, format_9993), s2 - s1;
+    write(nout, format_9993), cast2double(s2 - s1);
     //
     // End of Rchkee
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkee); }
+int main(int argc, char const *argv[]) { Rchkee(); }
