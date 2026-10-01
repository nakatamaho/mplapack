--- Ctsqr01.cpp_	2026-01-27 17:55:09.544577286 +0900
+++ Ctsqr01.cpp	2026-01-27 17:55:13.747659058 +0900
@@ -56,7 +56,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
@@ -157,7 +161,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", m, m, czero, one, r, m);
-        Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+        Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
         resid = Clansy("1", "Upper", m, r, m, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, m));
         //
@@ -288,7 +296,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", n, n, czero, one, lq, l);
-        Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), lq, l);
+        Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), lq, l);
         resid = Clansy("1", "Upper", n, lq, l, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, n));
         //
@@ -376,7 +376,7 @@
         // Compute |C*Q - C*Q| / |C|
         //
         Cgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
-        resid = Clange("1", n, m, df, n, rwork);
+        resid = Clange("1", m, n, cf, m, rwork);
         if (cnorm > zero) {
             result[5 - 1] = resid / (eps * max((INTEGER)1, n) * cnorm);
         } else {
