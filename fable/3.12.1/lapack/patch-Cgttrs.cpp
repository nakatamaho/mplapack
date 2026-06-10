--- a/mplapack/reference/Cgttrs.cpp
+++ b/mplapack/reference/Cgttrs.cpp
@@ -65,7 +65,7 @@
     INTEGER itrans = 0;
     if (notran) {
         itrans = 0;
-    } else if Mlsame(trans, "T") {
+    } else if (Mlsame(trans, "T")) {
         itrans = 1;
     } else {
         itrans = 2;
