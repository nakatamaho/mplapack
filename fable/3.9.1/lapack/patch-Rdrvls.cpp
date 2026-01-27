diff --git a/mplapack/test/lin/common/Rdrvls.cpp b/mplapack/test/lin/common/Rdrvls.cpp
index d48bef4c..f6bea5e8 100644
--- Rdrvls.cpp
+++ Rdrvls.cpp
@@ -94,8 +94,6 @@ void Rdrvls(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INT
     INTEGER nrows = 0;
     INTEGER ncols = 0;
     INTEGER ldwork = 0;
-    std::unique_ptr<REAL[]> work_storage(new REAL[lwork]);
-    REAL *work = work_storage.get();
     const REAL one = 1.0;
     const INTEGER ntests = 16;
     REAL result[ntests];
@@ -104,8 +102,6 @@ void Rdrvls(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INT
     INTEGER rank = 0;
     REAL normb = 0.0;
     INTEGER j = 0;
-    std::unique_ptr<INTEGER[]> iwork_storage(new INTEGER[liwork]);
-    INTEGER *iwork = iwork_storage.get();
     //
     static const char *format_9999 = "(' TRANS=''',a1,''', M=',i5,', N=',i5,', NRHS=',i4,', NB=',i4,', type',"
                                      "i2,', test(',i2,')=',g12.5)";
@@ -232,6 +228,10 @@ void Rdrvls(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INT
     //
     lwlsy = lwork;
     //
+    std::unique_ptr<REAL[]> work_storage(new REAL[lwork]);
+    REAL *work = work_storage.get();
+    std::unique_ptr<INTEGER[]> iwork_storage(new INTEGER[liwork]);
+    INTEGER *iwork = iwork_storage.get();
     for (im = 1; im <= nm; im = im + 1) {
         m = mval[im - 1];
         lda = max((INTEGER)1, m);
