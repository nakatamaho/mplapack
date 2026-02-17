--- Claror.cpp-	2026-01-18 08:42:53.359927619 +0900
+++ Claror.cpp	2026-01-18 08:42:59.389003982 +0900
@@ -42,6 +42,11 @@
 
 #include <mplapack_matgen.h>
 
+#if defined ___MPLAPACK_BUILD_WITH_DD___
+#pragma GCC push_options
+#pragma GCC optimize("O0")
+#endif
+
 void Claror(fem::str_cref side, fem::str_cref init, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER (&iseed)[4], COMPLEX *x, INTEGER &info) {
     //
     info = 0;
@@ -110,7 +115,7 @@
     COMPLEX csign = 0.0;
     COMPLEX xnorms = 0.0;
     REAL factor = 0.0;
-    const REAL toosml = 0.00000000000000000001;
+    const REAL toosml = 1.0e-20;
     const REAL one = 1.0;
     for (ixfrm = 2; ixfrm <= nxfrm; ixfrm = ixfrm + 1) {
         kbeg = nxfrm - ixfrm + 1;
@@ -140,7 +145,7 @@
         } else {
             factor = one / factor;
         }
-        x[kbeg - 1] += xnorms;
+        x[kbeg - 1] += xnorms; //this somehow doesn't work properly with GCC + libqd
         //
         // Apply Householder transformation to A
         //
@@ -202,3 +207,7 @@
     // End of Claror
     //
 }
+
+#if defined ___MPLAPACK_BUILD_WITH_DD___
+#pragma GCC pop_options
+#endif

@@ -106,7 +111,7 @@
     COMPLEX csign = 0.0;
     COMPLEX xnorms = 0.0;
     REAL factor = 0.0;
-    const REAL toosml = 0.00000000000000000001;
+    const REAL toosml = 1.0e-20;
     const REAL one = 1.0;
     for (ixfrm = 2; ixfrm <= nxfrm; ixfrm = ixfrm + 1) {
         kbeg = nxfrm - ixfrm + 1;
@@ -136,7 +141,7 @@
         } else {
             factor = one / factor;
         }
-        x[kbeg - 1] += xnorms;
+        x[kbeg - 1] += xnorms; //this somehow doesn't work properly with GCC + libqd
         //
         // Apply Householder transformation to A
         //
@@ -198,3 +203,7 @@
     // End of Claror
     //
 }
+
+#if defined ___MPLAPACK_BUILD_WITH_DD___
+#pragma GCC pop_options
+#endif

