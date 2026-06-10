--- Cdrvrf3.cpp.orig	2026-03-24 08:47:45.367979010 +0900
+++ Cdrvrf3.cpp	2026-03-24 08:56:12.069763248 +0900
@@ -96,6 +96,7 @@
     INTEGER j = 0;
     const INTEGER ntests = 1;
     REAL result[ntests];
+    const REAL two = 2.0;
     for (iim = 1; iim <= nn; iim = iim + 1) {
         //
         m = nval[iim - 1];
@@ -186,7 +187,7 @@
                                         if (Mlsame(diag.elems, "U")) {
                                             for (j = 1; j <= na; j = j + 1) {
                                                 for (i = 1; i <= j; i = i + 1) {
-                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] / (2.0 * a[(j - 1) + (j - 1) * lda]);
+                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] / (two * a[(j - 1) + (j - 1) * lda]);
                                                 }
                                             }
                                         }
@@ -206,7 +207,7 @@
                                         if (Mlsame(diag.elems, "U")) {
                                             for (i = 1; i <= na; i = i + 1) {
                                                 for (j = 1; j <= i; j = j + 1) {
-                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] / (2.0 * a[(i - 1) + (i - 1) * lda]);
+                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] / (two * a[(i - 1) + (i - 1) * lda]);
                                                 }
                                             }
                                         }
