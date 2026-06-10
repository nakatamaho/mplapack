--- Rchkec.cpp~	2026-01-22 11:26:12.550390264 +0900
+++ Rchkec.cpp	2026-01-22 11:27:03.798256608 +0900
@@ -181,7 +181,7 @@
     //
     REAL rtgexc = 0.0;
     INTEGER ltgexc = 0;
-    INTEGER ntgexc = 0;
+    INTEGER ntgexc[3]; //bug, fixed in lapack 3.12.1 correctly (dget40 is also correcte) to ninfo(2)
     INTEGER ktgexc = 0;
     Rget40(rtgexc, ltgexc, ntgexc, ktgexc, nin);
     if (rtgexc > thresh) {
