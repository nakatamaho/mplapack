--- Cqrt04.cpp_	2026-01-21 17:36:05.240536817 +0900
+++ Cqrt04.cpp	2026-01-21 17:36:08.823613729 +0900
@@ -94,7 +95,7 @@
     // Compute |I - Q'*Q| and store in RESULT(2)
     //
     Claset("Full", m, m, czero, one, r, m);
-    Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+    Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
