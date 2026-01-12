--- a/mplapack/reference/Ctrexc.cpp
+++ b/mplapack/reference/Ctrexc.cpp
@@ -91,7 +91,7 @@
     REAL cs = 0.0;
     COMPLEX sn = 0.0;
     COMPLEX temp = 0.0;
-    for (k = ifst + m1; k <= ilst + m2; k = k + m3) {
+    for (k = ifst + m1; m3 >= 0 ? k <= ilst + m2 : k >= ilst + m2; k = k + m3) {
         //
         // Interchange the k-th and (k+1)-th diagonal elements.
         //
