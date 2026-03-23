--- Cunhr_col02.cpp~	2026-03-23 17:38:25.096466598 +0900
+++ Cunhr_col02.cpp	2026-03-23 17:41:54.104928213 +0900
@@ -52,7 +52,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     //
     // Dynamically allocate local arrays
     //
@@ -191,7 +191,7 @@
     // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
     //
     Claset("Full", m, m, czero, cone, r, m);
-    Cherk("U", "C", m, m, -cone.real(), q, m, cone.real(), r, m);
+    Cherk("U", "C", m, m, -(cone).real(), q, m, cone.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
