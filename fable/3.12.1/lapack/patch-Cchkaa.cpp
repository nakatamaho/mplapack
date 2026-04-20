--- Cchkaa.cpp_org	2026-03-24 09:16:12.593866505 +0900
+++ Cchkaa.cpp	2026-03-24 09:29:05.932657463 +0900
@@ -44,38 +44,27 @@
 #include <mplapack_lin.h>
 #include <memory>
 
-void program_zchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Cchkaa(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
-    static fem::str<10> intstr = "0123456789";
-    static REAL threq = 2.0;
-    INTEGER allocatestatus = 0;
+    fem::str<10> intstr = "0123456789";
+    REAL threq = 2.0;
     const INTEGER nmax = 132;
-    const INTEGER kdmax = nmax + (nmax + 1) / 4;
-    std::unique_ptr<COMPLEX[]> a_storage;
-    COMPLEX *a = nullptr;
-    const INTEGER maxrhs = 16;
-    std::unique_ptr<COMPLEX[]> b_storage;
-    COMPLEX *b = nullptr;
-    std::unique_ptr<COMPLEX[]> work_storage;
-    COMPLEX *work = nullptr;
-    std::unique_ptr<COMPLEX[]> e_storage;
-    COMPLEX *e = nullptr;
-    std::unique_ptr<REAL[]> s_storage;
-    REAL *s = nullptr;
-    std::unique_ptr<REAL[]> rwork_storage;
-    REAL *rwork = nullptr;
     REAL s1 = 0.0;
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
+    const INTEGER maxrhs = 16;
     INTEGER mval[maxin];
     INTEGER i = 0;
     INTEGER nn = 0;
@@ -107,10 +96,21 @@
     fem::str<2> c2;
     INTEGER ntypes = 0;
     bool dotype[matmax];
+    const INTEGER kdmax = nmax + (nmax + 1) / 4;
+    auto a_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ((kdmax + 1) * nmax) * 7));
+    COMPLEX *a = a_storage.get();
+    auto b_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, (nmax * maxrhs) * 4));
+    COMPLEX *b = b_storage.get();
+    auto work_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * (nmax + maxrhs + 10)));
+    COMPLEX *work = work_storage.get();
+    auto rwork_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 150 * nmax + 2 * maxrhs));
+    REAL *rwork = rwork_storage.get();
     INTEGER iwork[25 * nmax];
+    REAL s[2 * nmax];
     INTEGER la = 0;
     INTEGER lafac = 0;
     INTEGER piv[nmax];
+    COMPLEX e[nmax];
     REAL s2 = 0.0;
     INTEGER ldaw = (kdmax + 1) * nmax;
     INTEGER ldb = nmax * maxrhs;
@@ -120,8 +120,9 @@
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
@@ -130,43 +131,6 @@
     static const char *format_9989 = "(/,1x,a3,' routines were not tested')";
     static const char *format_9988 = "(/,1x,a3,' driver routines were not tested')";
     //
-    allocatestatus = 0;
-    a_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ((kdmax + 1) * nmax) * 7));
-    a = a_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
-    b_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (nmax * maxrhs) * 4));
-    b = b_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
-    work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, nmax * (nmax + maxrhs + 10)));
-    work = work_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
-    e_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, nmax));
-    e = e_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
-    s_storage = std::make_unique<REAL[]>(max((INTEGER)1, (2 * nmax)));
-    s = s_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    allocatestatus = 0;
-    rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, (150 * nmax + 2 * maxrhs)));
-    rwork = rwork_storage.get();
-    if (allocatestatus != 0) {
-        FEM_STOP("*** Not enough memory ***");
-    }
-    //
     s1 = dsecnd();
     lda = nmax;
     fatal = false;
@@ -177,8 +141,8 @@
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
@@ -477,7 +441,7 @@
     //
     // Check first character for correct precision.
     //
-    if (!Mlsame(c1.elems, "Zomplex precision")) {
+    if (!Mlsame(c1.elems, "Zomplex precision") && !Mlsame(c1.elems, "C")) {
         write(nout, format_9990), path;
         //
     } else if (nmats <= 0) {
@@ -1015,6 +979,19 @@
             write(nout, format_9989), path;
         }
         //
+    } else if (Mlsamen(2, c2.elems, "QK")) {
+        //
+        // QK: truncated QR factorization with pivoting
+        //
+        ntypes = 19;
+        Alareq(path, nmats, dotype, ntypes, nin, nout);
+        //
+        if (tstchk) {
+            Cchkqp3rk(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &s[1 - 1], &b[(4 - 1) * ldb], work, rwork, iwork, nout);
+        } else {
+            write(nout, format_9989), path;
+        }
+        //
     } else if (Mlsamen(2, c2.elems, "LS")) {
         //
         // LS:  Least squares drivers
@@ -1133,10 +1110,10 @@
     cmn.io.close(nin);
     s2 = dsecnd();
     write(nout, format_9998);
-    write(nout, format_9997), s2 - s1;
+    write(nout, format_9997), cast2double(s2 - s1);
     //
     // End of Cchkaa
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkaa); }
+int main(int argc, char const *argv[]) { Cchkaa(); }
