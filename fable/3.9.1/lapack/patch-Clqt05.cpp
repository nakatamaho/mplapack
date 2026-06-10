--- Clqt05.cpp_	2026-01-21 17:43:44.795521221 +0900
+++ Clqt05.cpp	2026-01-21 17:44:26.348432698 +0900
@@ -112,7 +113,7 @@
     // Compute |I - Q*Q'| and store in RESULT(2)
     //
     Claset("Full", n2, n2, czero, one, r, n2);
-    Cherk("U", "N", n2, n2, dreal(-one), q, n2, dreal(one), r, n2);
+    Cherk("U", "N", n2, n2, (-one).real(), q, n2, one.real(), r, n2);
     resid = Clansy("1", "Upper", n2, r, n2, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, n2));
     //
