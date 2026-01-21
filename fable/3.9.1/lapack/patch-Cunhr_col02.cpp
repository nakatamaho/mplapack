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
