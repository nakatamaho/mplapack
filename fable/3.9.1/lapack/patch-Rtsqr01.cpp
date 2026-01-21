--- Rtsqr01.cpp_	2026-01-21 17:15:14.171909217 +0900
+++ Rtsqr01.cpp	2026-01-21 17:15:19.692982683 +0900
@@ -57,7 +57,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
