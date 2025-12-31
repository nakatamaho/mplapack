diff --git a/mplapack/reference/Clasyf_aa.cpp b/mplapack/reference/Clasyf_aa.cpp
index 4e83336c..999d848b 100644
--- a/mplapack/reference/Clasyf_aa.cpp
+++ b/mplapack/reference/Clasyf_aa.cpp
@@ -129,7 +129,7 @@ void Clasyf_aa(const char *uplo, INTEGER const j1, INTEGER const m, INTEGER cons
             //
             // Apply symmetric pivot
             //
-            if ((i2 != 2) && (piv != 0)) {
+            if ((i2 != 2) && (piv != zero)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
@@ -273,7 +273,7 @@ void Clasyf_aa(const char *uplo, INTEGER const j1, INTEGER const m, INTEGER cons
             //
             // Apply symmetric pivot
             //
-            if ((i2 != 2) && (piv != 0)) {
+            if ((i2 != 2) && (piv != zero)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
