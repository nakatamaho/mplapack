--- a/mplapack/reference/Clahef_aa.cpp
+++ b/mplapack/reference/Clahef_aa.cpp
@@ -131,7 +131,7 @@
             //
             // Apply hermitian pivot
             //
-            if ((i2 != 2) && (piv != 0)) {
+            if ((i2 != 2) && (piv != zero)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
@@ -279,7 +279,7 @@
             //
             // Apply hermitian pivot
             //
-            if ((i2 != 2) && (piv != 0)) {
+            if ((i2 != 2) && (piv != zero)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
