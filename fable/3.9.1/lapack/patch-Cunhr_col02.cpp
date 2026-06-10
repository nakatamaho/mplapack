--- Cunhr_col02.cpp_	2026-01-21 17:28:05.610488555 +0900
+++ Cunhr_col02.cpp	2026-01-21 17:28:19.487769970 +0900
@@ -53,7 +53,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     //
     // Dynamically allocate local arrays
     //
@@ -151,7 +158,7 @@
     // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
     //
     Claset("Full", m, m, czero, cone, r, m);
-    Cherk("U", "C", m, m, -cone, q, m, cone, r, m);
+    Cherk("U", "C", m, m, (-cone).real(), q, m, cone.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
