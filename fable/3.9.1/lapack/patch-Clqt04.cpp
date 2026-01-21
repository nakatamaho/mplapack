--- Clqt04.cpp_	2026-01-21 17:25:43.242648604 +0900
+++ Clqt04.cpp	2026-01-21 17:25:46.573713967 +0900
@@ -80,8 +81,8 @@
     // Compute |L - A*Q'| / |A| and store in RESULT(1)
     //
     Cgemm("N", "C", m, n, n, -one, a, m, q, n, one, l, ll);
-    std::unique_ptr<COMPLEX[]> __rwork_storage(new COMPLEX[ll]);
-    COMPLEX *rwork = __rwork_storage.get();
+    std::unique_ptr<REAL[]> __rwork_storage(new REAL[ll]);
+    REAL *rwork = __rwork_storage.get();
     REAL anorm = Clange("1", m, n, a, m, rwork);
     REAL resid = Clange("1", m, n, l, ll, rwork);
     const REAL zero = 0.0;
@@ -94,7 +95,7 @@
     // Compute |I - Q'*Q| and store in RESULT(2)
     //
     Claset("Full", n, n, czero, one, l, ll);
-    Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), l, ll);
+    Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), l, ll);
     resid = Clansy("1", "Upper", n, l, ll, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, n));
     //
