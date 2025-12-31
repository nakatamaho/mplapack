diff --git a/mplapack/reference/iMieeeck.cpp b/mplapack/reference/iMieeeck.cpp
index 125347e3..3e59de82 100644
--- a/mplapack/reference/iMieeeck.cpp
+++ b/mplapack/reference/iMieeeck.cpp
@@ -36,20 +36,9 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-INTEGER iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
+INTEGER
+iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
     INTEGER return_value = 0;
-#if defined ___MPLAPACK_BUILD_WITH_GMP___
-    // GMP is not a natural extention to IEEE 754.
-    return 0;
-#endif
-#if defined ___MPLAPACK_BUILD_WITH_DD___
-    // DD does not comply IEEE 754.
-    return 0;
-#endif
-#if defined ___MPLAPACK_BUILD_WITH_QD___
-    // DD does not comply IEEE 754.
-    return 0;
-#endif
     return_value = 1;
     //
     REAL posinf = one / zero;
