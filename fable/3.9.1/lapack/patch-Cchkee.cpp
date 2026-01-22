--- Cchkee.cpp
+++ Cchkee.cpp
@@ -43,19 +43,15 @@ using fem::common;
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
 
-void program_zchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
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
@@ -87,9 +83,12 @@ void program_zchkee(int argc, char const *argv[]) {
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
@@ -160,13 +159,21 @@ void program_zchkee(int argc, char const *argv[]) {
     //
     //
     //
-    a = 0.0;
-    b = 0.0;
-    c = 0.0;
-    dc = 0.0;
+    constexpr std::size_t NA = (std::size_t)nmax * nmax * need;
+    constexpr std::size_t NB = (std::size_t)nmax * nmax * 5;
+    constexpr std::size_t NC = (std::size_t)ncmax * ncmax * ncmax * ncmax;
+    constexpr std::size_t ND = (std::size_t)nmax * 6;
+    COMPLEX a[NA];
+    COMPLEX b[NB];
+    COMPLEX c[NC];
+    COMPLEX dc[ND];
+    const COMPLEX zero = 0.0;
+    std::fill_n(a, NA, zero);
+    std::fill_n(b, NB, zero);
+    std::fill_n(c, NC, zero);
+    std::fill_n(dc, ND, zero);
     s1 = dsecnd();
     fatal = false;
-    nunit = nout;
 //
 // Return to here to read multiple sets of data
 //
@@ -296,8 +303,10 @@ statement_10:
         write(nout, format_9992), path;
         goto statement_380;
     }
-    ilaver(vers_major, vers_minor, vers_patch);
-    write(nout, "(/,' LAPACK VERSION ',i1,'.',i1,'.',i1)"), vers_major, vers_minor, vers_patch;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
+                "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
+        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     write(nout, "(/,' The following parameter values will be used:')");
     //
     // Read the number of values of M, P, and N.
@@ -1027,7 +1036,7 @@ statement_190:
             }
             write(nout, "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', INMIN=',i4,"
                         "', INWIN =',i4,', INIBL =',i4,', ISHFTS =',i4,', IACC22 =',i4)"),
-                c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], fem::max(11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+                c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
             Cchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &dc[(3 - 1) * nmax], work, lwork, rwork, iwork, logwrk, result, info);
             if (info != 0) {
                 write(nout, format_9980), "ZCHKHS", info;
@@ -1565,4 +1574,4 @@ statement_380:
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkee); }
+int main(int argc, char const *argv[]) { Cchkee(); }
