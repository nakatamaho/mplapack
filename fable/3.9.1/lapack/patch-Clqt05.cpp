--- Clqt05.cpp_	2026-01-21 17:43:44.795521221 +0900
+++ Clqt05.cpp	2026-01-21 17:44:26.348432698 +0900
@@ -98,8 +99,8 @@
     // Compute |L - A*Q*C| / |A| and store in RESULT(1)
     //
     Cgemm("N", "C", m, n2, n2, -one, a, m, q, n2, one, r, n2);
-    std::unique_ptr<COMPLEX[]> rwork_storage(new COMPLEX[n2]);
-    COMPLEX *rwork = rwork_storage.get();
+    std::unique_ptr<REAL[]> rwork_storage(new REAL[n2]);
+    REAL *rwork = rwork_storage.get();
     REAL anorm = Clange("1", m, n2, a, m, rwork);
     REAL resid = Clange("1", m, n2, r, n2, rwork);
     const REAL zero = 0.0;
@@ -112,7 +113,7 @@
     // Compute |I - Q*Q'| and store in RESULT(2)
     //
     Claset("Full", n2, n2, czero, one, r, n2);
-    Cherk("U", "N", n2, n2, dreal(-one), q, n2, dreal(one), r, n2);
+    Cherk("U", "N", n2, n2, (-one).real(), q, n2, one.real(), r, n2);
     resid = Clansy("1", "Upper", n2, r, n2, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, n2));
     //
