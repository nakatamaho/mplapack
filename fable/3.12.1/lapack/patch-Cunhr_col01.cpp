--- Cunhr_col01.cpp~	2026-03-23 17:38:25.450474154 +0900
+++ Cunhr_col01.cpp	2026-03-23 17:40:36.440270261 +0900
@@ -52,7 +52,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     //
     // Dynamically allocate local arrays
     //
@@ -230,7 +230,7 @@
     // Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)
     //
     Claset("Full", m, m, czero, cone, r, m);
-    Cherk("U", "C", m, m, -cone.real(), q, m, cone.real(), r, m);
+    Cherk("U", "C", m, m, -(cone).real(), q, m, cone.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
