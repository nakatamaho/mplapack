diff --git a/mplapack/reference/Cgttrs.cpp b/mplapack/reference/Cgttrs.cpp
index 78edcf62..adbb6f3b 100644
--- a/mplapack/reference/Cgttrs.cpp
+++ b/mplapack/reference/Cgttrs.cpp
@@ -65,7 +65,7 @@ void Cgttrs(const char *trans, INTEGER const n, INTEGER const nrhs, COMPLEX *dl,
     INTEGER itrans = 0;
     if (notran) {
         itrans = 0;
-    } else if Mlsame(trans, "T") {
+    } else if (Mlsame(trans, "T")) {
         itrans = 1;
     } else {
         itrans = 2;
