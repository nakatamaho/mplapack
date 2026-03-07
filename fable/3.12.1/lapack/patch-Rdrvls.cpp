--- Rdrvls.cpp_	2026-01-27 17:39:46.555713071 +0900
+++ Rdrvls.cpp	2026-01-27 17:41:17.547701216 +0900
@@ -95,7 +95,7 @@
     INTEGER ncols = 0;
     INTEGER ldwork = 0;
     std::unique_ptr<REAL[]> work_storage(new REAL[lwork]);
-    REAL *work = work_storage.get();
+    REAL *work = nullptr;
     const REAL one = 1.0;
     const INTEGER ntests = 16;
     REAL result[ntests];
@@ -105,7 +105,7 @@
     REAL normb = 0.0;
     INTEGER j = 0;
     std::unique_ptr<INTEGER[]> iwork_storage(new INTEGER[liwork]);
-    INTEGER *iwork = iwork_storage.get();
+    INTEGER *iwork = nullptr;
     //
     static const char *format_9999 = "(' TRANS=''',a1,''', M=',i5,', N=',i5,', NRHS=',i4,', NB=',i4,', type',"
                                      "i2,', test(',i2,')=',g12.5)";
@@ -232,6 +232,10 @@
     //
     lwlsy = lwork;
     //
+    work_storage.reset(new REAL[max((INTEGER)1, lwork)]);
+    work = work_storage.get();
+    iwork_storage.reset(new INTEGER[max((INTEGER)1, liwork)]);
+    iwork = iwork_storage.get();
     for (im = 1; im <= nm; im = im + 1) {
         m = mval[im - 1];
         lda = max((INTEGER)1, m);
