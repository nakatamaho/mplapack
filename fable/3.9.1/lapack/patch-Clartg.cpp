diff --git a/mplapack/reference/Clartg.cpp b/mplapack/reference/Clartg.cpp
index 239bcd84..e23392bd 100644
--- a/mplapack/reference/Clartg.cpp
+++ b/mplapack/reference/Clartg.cpp
@@ -36,6 +36,12 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL abssq(COMPLEX ff) {
+    REAL temp;
+    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
+    return temp;
+}
+
 void Clartg(COMPLEX const f, COMPLEX const g, REAL &cs, COMPLEX &sn, COMPLEX &r) {
     COMPLEX ff = 0.0;
     REAL safmin = 0.0;
