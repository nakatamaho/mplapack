--- a/mplapack/reference/Clasyf_aa.cpp
+++ b/mplapack/reference/Clasyf_aa.cpp
@@ -129,7 +129,7 @@
             //
             // Apply symmetric pivot
             //
-            if ((i2 != 2) && (piv != 0)) {
+            if ((i2 != 2) && (piv != zero)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
@@ -273,7 +273,7 @@
             //
             // Apply symmetric pivot
             //
-            if ((i2 != 2) && (piv != 0)) {
+            if ((i2 != 2) && (piv != zero)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
