--- Cqrt05.cpp_	2026-01-21 17:23:09.003691019 +0900
+++ Cqrt05.cpp	2026-01-21 17:23:20.006897091 +0900
@@ -112,7 +112,7 @@
     // Compute |I - Q'*Q| and store in RESULT(2)
     //
     Claset("Full", m2, m2, czero, one, r, m2);
-    Cherk("U", "C", m2, m2, dreal(-one), q, m2, dreal(one), r, m2);
+    Cherk("U", "C", m2, m2, -one.real(), q, m2, one.real(), r, m2);
     resid = Clansy("1", "Upper", m2, r, m2, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m2));
     //
