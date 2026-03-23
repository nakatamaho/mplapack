diff --git a/mplapack/test/lin/common/Cchkaa.cpp b/mplapack/test/lin/common/Cchkaa.cpp
index 642064ad..1fc0229e 100644
--- a/mplapack/test/lin/common/Cchkaa.cpp
+++ b/mplapack/test/lin/common/Cchkaa.cpp
@@ -50,8 +50,23 @@ void Cchkaa(void) {
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
     static REAL threq = 2.0;
-    REAL s1 = 0.0;
+    INTEGER allocatestatus = 0;
     const INTEGER nmax = 132;
+    const INTEGER kdmax = nmax + (nmax + 1) / 4;
+    std::unique_ptr<COMPLEX[]> a_storage;
+    COMPLEX *a = nullptr;
+    const INTEGER maxrhs = 16;
+    std::unique_ptr<COMPLEX[]> b_storage;
+    COMPLEX *b = nullptr;
+    std::unique_ptr<COMPLEX[]> work_storage;
+    COMPLEX *work = nullptr;
+    std::unique_ptr<COMPLEX[]> e_storage;
+    COMPLEX *e = nullptr;
+    std::unique_ptr<REAL[]> s_storage;
+    REAL *s = nullptr;
+    std::unique_ptr<REAL[]> rwork_storage;
+    REAL *rwork = nullptr;
+    REAL s1 = 0.0;
     INTEGER lda = 0;
     bool fatal = false;
     const INTEGER nin = 5;
@@ -70,7 +85,6 @@ void Cchkaa(void) {
     INTEGER nval[maxin];
     INTEGER nns = 0;
     INTEGER nsval[maxin];
-    const INTEGER maxrhs = 16;
     INTEGER nnb = 0;
     INTEGER nbval[maxin];
     INTEGER nnb2 = 0;
@@ -96,21 +110,10 @@ void Cchkaa(void) {
     fem::str<2> c2;
     INTEGER ntypes = 0;
     bool dotype[matmax];
-    const INTEGER kdmax = nmax + (nmax + 1) / 4;
-    auto a_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ((kdmax + 1) * nmax) * 7));
-    COMPLEX *a = a_storage.get();
-    auto b_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, (nmax * maxrhs) * 4));
-    COMPLEX *b = b_storage.get();
-    auto work_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * (nmax + maxrhs + 10)));
-    COMPLEX *work = work_storage.get();
-    auto rwork_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 150 * nmax + 2 * maxrhs));
-    REAL *rwork = rwork_storage.get();
     INTEGER iwork[25 * nmax];
-    REAL s[2 * nmax];
     INTEGER la = 0;
     INTEGER lafac = 0;
     INTEGER piv[nmax];
-    COMPLEX e[nmax];
     REAL s2 = 0.0;
     INTEGER ldaw = (kdmax + 1) * nmax;
     INTEGER ldb = nmax * maxrhs;
@@ -131,7 +134,6 @@ void Cchkaa(void) {
     static const char *format_9989 = "(/,1x,a3,' routines were not tested')";
     static const char *format_9988 = "(/,1x,a3,' driver routines were not tested')";
     //
-    //
     s1 = dsecnd();
     lda = nmax;
     fatal = false;
@@ -967,6 +969,32 @@ statement_130:
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
