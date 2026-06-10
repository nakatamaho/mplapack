--- Cdrvls.cpp_	2026-01-27 17:40:03.077075407 +0900
+++ Cdrvls.cpp	2026-01-27 17:40:15.003336573 +0900
@@ -100,23 +100,23 @@
     INTEGER ncols = 0;
     INTEGER ldwork = 0;
     std::unique_ptr<COMPLEX[]> work_storage(new COMPLEX[lwork]);
-    COMPLEX *work = work_storage.get();
+    COMPLEX *work = nullptr;
     const REAL one = 1.0;
     const COMPLEX cone = COMPLEX(1.0, 0.0);
     const COMPLEX czero = COMPLEX(0.0, 0.0);
     std::unique_ptr<REAL[]> rwork_storage(new REAL[lrwork]);
-    REAL *rwork = rwork_storage.get();
+    REAL *rwork = nullptr;
     const INTEGER ntests = 16;
     REAL result[ntests];
     INTEGER k = 0;
     INTEGER imb = 0;
     std::unique_ptr<REAL[]> work2_storage(new REAL[2 * lwork]);
-    REAL *work2 = work2_storage.get();
+    REAL *work2 = nullptr;
     INTEGER rank = 0;
     REAL normb = 0.0;
     INTEGER j = 0;
     std::unique_ptr<INTEGER[]> iwork_storage(new INTEGER[liwork]);
-    INTEGER *iwork = iwork_storage.get();
+    INTEGER *iwork = nullptr;
     //
     static const char *format_9999 = "(' TRANS=''',a1,''', M=',i5,', N=',i5,', NRHS=',i4,', NB=',i4,', type',"
                                      "i2,', test(',i2,')=',g12.5)";
@@ -245,6 +245,14 @@
     }
     //
     lwlsy = lwork;
+    work_storage.reset(new COMPLEX[max((INTEGER)1, lwork)]);
+    work = work_storage.get();
+    work2_storage.reset(new REAL[max((INTEGER)1, (INTEGER)2 * lwork)]);
+    work2 = work2_storage.get();
+    iwork_storage.reset(new INTEGER[max((INTEGER)1, liwork)]);
+    iwork = iwork_storage.get();
+    rwork_storage.reset(new REAL[max((INTEGER)1, lrwork)]);
+    rwork = rwork_storage.get();
     //
     for (im = 1; im <= nm; im = im + 1) {
         m = mval[im - 1];
