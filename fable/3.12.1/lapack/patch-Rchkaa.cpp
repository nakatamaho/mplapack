diff --git a/mplapack/test/lin/common/Rchkaa.cpp b/mplapack/test/lin/common/Rchkaa.cpp
index 8cb12ec6..e83607cb 100644
--- a/mplapack/test/lin/common/Rchkaa.cpp
+++ b/mplapack/test/lin/common/Rchkaa.cpp
@@ -51,8 +51,23 @@ void Rchkaa(void) {
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
     static REAL threq = 2.0;
-    REAL s1 = 0.0;
+    INTEGER allocatestatus = 0;
     const INTEGER nmax = 132;
+    const INTEGER kdmax = nmax + (nmax + 1) / 4;
+    std::unique_ptr<REAL[]> a_storage;
+    REAL *a = nullptr;
+    const INTEGER maxrhs = 16;
+    std::unique_ptr<REAL[]> b_storage;
+    REAL *b = nullptr;
+    std::unique_ptr<REAL[]> work_storage;
+    REAL *work = nullptr;
+    std::unique_ptr<REAL[]> e_storage;
+    REAL *e = nullptr;
+    std::unique_ptr<REAL[]> s_storage;
+    REAL *s = nullptr;
+    std::unique_ptr<REAL[]> rwork_storage;
+    REAL *rwork = nullptr;
+    REAL s1 = 0.0;
     INTEGER lda = 0;
     bool fatal = false;
     const INTEGER nin = 5;
@@ -71,7 +86,6 @@ void Rchkaa(void) {
     INTEGER nval[maxin];
     INTEGER nns = 0;
     INTEGER nsval[maxin];
-    const INTEGER maxrhs = 16;
     INTEGER nnb = 0;
     INTEGER nbval[maxin];
     INTEGER nnb2 = 0;
@@ -97,21 +111,10 @@ void Rchkaa(void) {
     INTEGER nrhs = 0;
     INTEGER ntypes = 0;
     bool dotype[matmax];
-    const INTEGER kdmax = nmax + (nmax + 1) / 4;
-    auto a_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, ((kdmax + 1) * nmax) * 7));
-    REAL *a = a_storage.get();
-    auto b_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * maxrhs) * 4));
-    REAL *b = b_storage.get();
-    auto work_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * (3 * nmax + maxrhs + 30)));
-    REAL *work = work_storage.get();
-    auto rwork_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 5 * nmax + 2 * maxrhs));
-    REAL *rwork = rwork_storage.get();
     INTEGER iwork[25 * nmax];
-    REAL s[2 * nmax];
     INTEGER la = 0;
     INTEGER lafac = 0;
     INTEGER piv[nmax];
-    REAL e[nmax];
     REAL s2 = 0.0;
     INTEGER ldaw = (kdmax + 1) * nmax;
     INTEGER ldb = nmax * maxrhs;
@@ -132,7 +135,6 @@ void Rchkaa(void) {
     static const char *format_9989 = "(/,1x,a3,' routines were not tested')";
     static const char *format_9988 = "(/,1x,a3,' driver routines were not tested')";
     //
-    //
     s1 = dsecnd();
     lda = nmax;
     fatal = false;
@@ -827,6 +829,32 @@ statement_130:
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
+            Rchkqp3rk(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, iwork, nout);
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
+            Rchkqp3rk(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, iwork, nout);
+        } else {
+            write(nout, format_9989), path;
+        }
+        //
     } else if (Mlsamen(2, c2.elems, "TZ")) {
         //
         // TZ:  Trapezoidal matrix
