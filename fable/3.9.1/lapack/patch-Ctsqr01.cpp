--- Ctsqr01.cpp_	2026-01-21 17:33:58.844837374 +0900
+++ Ctsqr01.cpp	2026-01-21 17:34:07.203014994 +0900
@@ -57,7 +59,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
@@ -100,8 +102,8 @@
     COMPLEX *q = __q_storage.get();
     std::unique_ptr<COMPLEX[]> __r_storage(new COMPLEX[m * l]);
     COMPLEX *r = __r_storage.get();
-    std::unique_ptr<COMPLEX[]> __rwork_storage(new COMPLEX[l]);
-    COMPLEX *rwork = __rwork_storage.get();
+    std::unique_ptr<REAL[]> __rwork_storage(new REAL[l]);
+    REAL *rwork = __rwork_storage.get();
     REAL anorm = 0.0;
     REAL resid = 0.0;
     const REAL zero = 0.0;
@@ -158,7 +161,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", m, m, czero, one, r, m);
-        Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+        Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
         resid = Clansy("1", "Upper", m, r, m, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, m));
         //
@@ -289,7 +293,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", n, n, czero, one, lq, l);
-        Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), lq, l);
+        Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), lq, l);
         resid = Clansy("1", "Upper", n, lq, l, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, n));
         //
