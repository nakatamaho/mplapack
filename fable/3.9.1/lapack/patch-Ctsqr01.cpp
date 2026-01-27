--- Ctsqr01.cpp_	2026-01-27 17:55:09.544577286 +0900
+++ Ctsqr01.cpp	2026-01-27 17:55:13.747659058 +0900
@@ -56,7 +56,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
@@ -90,9 +90,9 @@
     std::unique_ptr<COMPLEX[]> df_storage(new COMPLEX[n * m]);
     COMPLEX *df = df_storage.get();
     std::unique_ptr<COMPLEX[]> t_storage(new COMPLEX[tsize]);
-    COMPLEX *t = t_storage.get();
+    COMPLEX *t = nullptr;
     std::unique_ptr<COMPLEX[]> work_storage(new COMPLEX[lwork]);
-    COMPLEX *work = work_storage.get();
+    COMPLEX *work = nullptr;
     const COMPLEX czero = COMPLEX(0.0, 0.0);
     const COMPLEX one = COMPLEX(1.0, 0.0);
     std::unique_ptr<COMPLEX[]> q_storage(new COMPLEX[l * l]);
@@ -129,6 +129,10 @@
         lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
         Cgemqr("R", "C", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
         lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
+        t_storage.reset(new COMPLEX[max((INTEGER)1, tsize)]);
+        t = t_storage.get();
+        work_storage.reset(new COMPLEX[max((INTEGER)1, lwork)]);
+        work = work_storage.get();
         srnamt = "Cgeqr";
         Cgeqr(m, n, af, m, t, tsize, work, lwork, info);
         //
@@ -157,7 +161,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", m, m, czero, one, r, m);
-        Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+        Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
         resid = Clansy("1", "Upper", m, r, m, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, m));
         //
@@ -260,6 +264,10 @@
         lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
         Cgemlq("R", "C", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
         lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
+        t_storage.reset(new COMPLEX[max((INTEGER)1, tsize)]);
+        t = t_storage.get();
+        work_storage.reset(new COMPLEX[max((INTEGER)1, lwork)]);
+        work = work_storage.get();
         srnamt = "Cgelq";
         Cgelq(m, n, af, m, t, tsize, work, lwork, info);
         //
@@ -288,7 +296,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", n, n, czero, one, lq, l);
-        Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), lq, l);
+        Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), lq, l);
         resid = Clansy("1", "Upper", n, lq, l, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, n));
         //
