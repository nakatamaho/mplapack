--- Rtsqr01.cpp_	2026-01-21 20:55:15.937518627 +0900
+++ Rtsqr01.cpp	2026-01-21 20:55:15.943518710 +0900
@@ -45,6 +45,7 @@
 #include <memory>
 
 void Rtsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *result) {
+    common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
     // TEST TALL SKINNY OR SHORT WIDE
     //
@@ -56,7 +57,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
