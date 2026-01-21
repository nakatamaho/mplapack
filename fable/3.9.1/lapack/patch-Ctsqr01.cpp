--- Ctsqr01.cpp_	2026-01-21 20:55:15.883517869 +0900
+++ Ctsqr01.cpp	2026-01-21 20:55:15.889517953 +0900
@@ -45,6 +45,7 @@
 #include <memory>
 
 void Ctsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *result) {
+    common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
     // TEST TALL SKINNY OR SHORT WIDE
     //
@@ -56,7 +57,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
@@ -99,8 +100,8 @@
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
@@ -157,7 +158,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", m, m, czero, one, r, m);
-        Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+        Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
         resid = Clansy("1", "Upper", m, r, m, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, m));
         //
@@ -288,7 +289,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", n, n, czero, one, lq, l);
-        Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), lq, l);
+        Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), lq, l);
         resid = Clansy("1", "Upper", n, lq, l, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, n));
         //
