--- Cchkgg.cpp~	2026-01-22 16:10:28.243935401 +0900
+++ Cchkgg.cpp	2026-01-22 16:13:10.532251318 +0900
@@ -115,7 +115,7 @@
         }
     }
     //
-    lwkopt = max(2 * nmax * nmax, 4 * nmax, 1);
+    lwkopt = max(2 * nmax * nmax, 4 * nmax, (INTEGER)1);
     //
     // Check for errors
     //
