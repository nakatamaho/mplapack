--- Rorhr_col01.cpp~	2026-01-21 17:14:19.743214542 +0900
+++ Rorhr_col01.cpp	2026-01-21 17:17:18.412681651 +0900
@@ -53,7 +53,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     //
     // Dynamically allocate local arrays
     //
