--- Clatm5.cpp_	2026-01-17 18:51:56.074139937 +0900
+++ Clatm5.cpp	2026-01-17 18:52:06.957536960 +0900
@@ -47,6 +47,7 @@
     const COMPLEX half = COMPLEX(0.5, 0.0);
     const COMPLEX twenty = COMPLEX(20.0, 0.0);
     const COMPLEX two = COMPLEX(2.0, 0.0);
+    const REAL rtwo = 2.0;
     INTEGER k = 0;
     COMPLEX reeps = 0.0;
     COMPLEX imeps = 0.0;
@@ -199,9 +200,9 @@
             } else {
                 a[(i - 1) + (i - 1) * lda] = one;
                 if (mod(i, 2) != 0 && i < m) {
-                    a[(i - 1) + ((i + 1) - 1) * lda] = imeps * 2;
+                    a[(i - 1) + ((i + 1) - 1) * lda] = imeps * rtwo;
                 } else if (i > 1) {
-                    a[(i - 1) + ((i - 1) - 1) * lda] = -imeps * 2;
+                    a[(i - 1) + ((i - 1) - 1) * lda] = -imeps * rtwo;
                 }
             }
         }
@@ -232,9 +233,9 @@
             } else {
                 b[(i - 1) + (i - 1) * ldb] = one - reeps;
                 if (mod(i, 2) != 0 && i < n) {
-                    b[(i - 1) + ((i + 1) - 1) * ldb] = imeps * 2;
+                    b[(i - 1) + ((i + 1) - 1) * ldb] = imeps * rtwo;
                 } else if (i > 1) {
-                    b[(i - 1) + ((i - 1) - 1) * ldb] = -imeps * 2;
+                    b[(i - 1) + ((i - 1) - 1) * ldb] = -imeps * rtwo;
                 }
             }
         }
