--- Rtsqr01.cpp_	2026-01-21 20:55:15.937518627 +0900
+++ Rtsqr01.cpp	2026-01-21 20:55:15.943518710 +0900
@@ -56,7 +57,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
@@ -90,9 +90,9 @@ void Rtsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const
     std::unique_ptr<REAL[]> df_storage(new REAL[n * m]);
     REAL *df = df_storage.get();
     std::unique_ptr<REAL[]> t_storage(new REAL[tsize]);
-    REAL *t = t_storage.get();
+    REAL *t = nullptr;
     std::unique_ptr<REAL[]> work_storage(new REAL[lwork]);
-    REAL *work = work_storage.get();
+    REAL *work = nullptr;
     const REAL zero = 0.0;
     const REAL one = 1.0;
     std::unique_ptr<REAL[]> q_storage(new REAL[l * l]);
@@ -128,6 +128,10 @@ void Rtsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const
         lwork = max(lwork, castINTEGER(workquery[1 - 1]));
         Rgemqr("R", "T", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
         lwork = max(lwork, castINTEGER(workquery[1 - 1]));
+        t_storage.reset(new REAL[max((INTEGER)1, tsize)]);
+        t = t_storage.get();
+        work_storage.reset(new REAL[max((INTEGER)1, lwork)]);
+        work = work_storage.get();
         srnamt = "Rgeqr";
         Rgeqr(m, n, af, m, t, tsize, work, lwork, info);
         //
@@ -259,6 +263,10 @@ void Rtsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const
         lwork = max(lwork, castINTEGER(workquery[1 - 1]));
         Rgemlq("R", "T", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
         lwork = max(lwork, castINTEGER(workquery[1 - 1]));
+        t_storage.reset(new REAL[max((INTEGER)1, tsize)]);
+        t = t_storage.get();
+        work_storage.reset(new REAL[max((INTEGER)1, lwork)]);
+        work = work_storage.get();
         srnamt = "Rgelq";
         Rgelq(m, n, af, m, t, tsize, work, lwork, info);
         //
