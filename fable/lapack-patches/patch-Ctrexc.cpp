diff --git a/mplapack/reference/Ctrexc.cpp b/mplapack/reference/Ctrexc.cpp
index e7445e1a..9242e20b 100644
--- a/mplapack/reference/Ctrexc.cpp
+++ b/mplapack/reference/Ctrexc.cpp
@@ -84,7 +84,7 @@ void Ctrexc(const char *compq, INTEGER const n, COMPLEX *t, INTEGER const ldt, C
     REAL cs = 0.0;
     COMPLEX sn = 0.0;
     COMPLEX temp = 0.0;
-    for (k = ifst + m1; m3 >= 0 ? k <= ilst + m2 : k >= ilst + m2; k = k + m3) {
+    for (k = ifst + m1; k <= ilst + m2; k = k + m3) {
         //
         // Interchange the k-th and (k+1)-th diagonal elements.
         //
