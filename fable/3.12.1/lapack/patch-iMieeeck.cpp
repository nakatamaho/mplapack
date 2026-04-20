--- a/mplapack/reference/iMieeeck.cpp
+++ b/mplapack/reference/iMieeeck.cpp
@@ -36,9 +36,27 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-INTEGER
-iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
+INTEGER iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
     INTEGER return_value = 0;
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+    // GMP uses arbitrary-precision integers/rationals internally and does not
+    // implement IEEE 754 semantics: division by zero does not yield +Inf/-Inf,
+    // and NaN/signed-zero behavior is undefined.  The runtime checks below
+    // would invoke undefined behavior on GMP arithmetic, so return 0 here.
+    return 0;
+#endif
+#if defined ___MPLAPACK_BUILD_WITH_DD___
+    // DD (double-double) arithmetic does not comply with IEEE 754: it lacks
+    // proper handling of infinities, NaN, and signed zero.  The runtime checks
+    // below would invoke undefined behavior on DD arithmetic, so return 0 here.
+    return 0;
+#endif
+#if defined ___MPLAPACK_BUILD_WITH_QD___
+    // QD (quad-double) arithmetic does not comply with IEEE 754: it lacks
+    // proper handling of infinities, NaN, and signed zero.  The runtime checks
+    // below would invoke undefined behavior on QD arithmetic, so return 0 here.
+    return 0;
+#endif
     return_value = 1;
     //
     REAL posinf = one / zero;
