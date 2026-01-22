--- Rchkee.cpp
+++ Rchkee.cpp
@@ -42,20 +42,17 @@ using fem::common;
 
 #include <mplapack_matgen.h>
 #include <mplapack_eig.h>
+#include <mplapack_debug.h>
 
-void program_dchkee(int argc, char const *argv[]) {
-    common cmn(argc, argv);
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
@@ -87,9 +84,12 @@ void program_dchkee(int argc, char const *argv[]) {
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
@@ -156,13 +156,21 @@ void program_dchkee(int argc, char const *argv[]) {
     //
     //
     //
-    a = 0.0;
-    b = 0.0;
-    c = 0.0;
-    d = 0.0;
+    constexpr std::size_t NA = (std::size_t)nmax * nmax * need;
+    constexpr std::size_t NB = (std::size_t)nmax * nmax * 5;
+    constexpr std::size_t NC = (std::size_t)ncmax * ncmax * ncmax * ncmax;
+    constexpr std::size_t ND = (std::size_t)nmax * 12;
+    REAL a[nmax * nmax * need];
+    REAL b[nmax * nmax * 5];
+    REAL c[ncmax * ncmax * ncmax * ncmax];
+    REAL d[nmax * 12];
+    const REAL zero = 0.0;
+    std::fill_n(a, NA, zero);
+    std::fill_n(b, NB, zero);
+    std::fill_n(c, NC, zero);
+    std::fill_n(d, ND, zero);
     s1 = dsecnd();
     fatal = false;
-    nunit = nout;
 //
 // Return to here to read multiple sets of data
 //
@@ -296,8 +304,10 @@ statement_10:
         write(nout, format_9992), path;
         goto statement_10;
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
@@ -1028,7 +1038,7 @@ statement_190:
             }
             write(nout, "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', INMIN=',i4,"
                         "', INWIN =',i4,', INIBL =',i4,', ISHFTS =',i4,', IACC22 =',i4)"),
-                c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], fem::max(11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
+                c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
             Rchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &d[0], &d[(2 - 1) * nmax], &d[(3 - 1) * nmax], &d[(4 - 1) * nmax], &d[(5 - 1) * nmax], &d[(6 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &d[(7 - 1) * nmax], work, lwork, iwork, logwrk, result, info);
             if (info != 0) {
                 write(nout, format_9980), "DCHKHS", info;
@@ -1563,4 +1573,4 @@ statement_380:
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkee); }
+int main(int argc, char const *argv[]) { Rchkee(); }
