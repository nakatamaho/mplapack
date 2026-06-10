--- Rchkgg.cpp~	2026-01-22 16:10:26.564900901 +0900
+++ Rchkgg.cpp	2026-01-22 16:12:47.646785759 +0900
@@ -114,7 +114,7 @@
     // Maximum blocksize and shift -- we assume that blocksize and number
     // of shifts are monotone increasing functions of N.
     //
-    lwkopt = max(6 * nmax, 2 * nmax * nmax, 1);
+    lwkopt = max(6 * nmax, 2 * nmax * nmax, (INTEGER)1);
     //
     // Check for errors
     //
