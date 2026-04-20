--- Rlaror.cpp-	2026-01-18 09:20:01.266785753 +0900
+++ Rlaror.cpp	2026-01-18 09:20:07.705936558 +0900
@@ -100,7 +100,7 @@
     REAL xnorm = 0.0;
     REAL xnorms = 0.0;
     REAL factor = 0.0;
-    const REAL toosml = 0.00000000000000000001;
+    const REAL toosml = 1.0e-20;
     for (ixfrm = 2; ixfrm <= nxfrm; ixfrm = ixfrm + 1) {
         kbeg = nxfrm - ixfrm + 1;
         //
