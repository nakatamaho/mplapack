--- a/mplapack/reference/iMieeeck.cpp
+++ b/mplapack/reference/iMieeeck.cpp
@@ -36,19 +36,37 @@
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
-    if (posinf <= one) {
+    if (!Risinf(posinf) || posinf <= zero) {
         return_value = 0;
         return return_value;
     }
     //
     REAL neginf = -one / zero;
-    if (neginf >= zero) {
+    if (!Risinf(neginf) || neginf >= zero) {
         return_value = 0;
         return return_value;
     }
@@ -61,7 +79,7 @@ INTEGER iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
     }
     //
     neginf = one / negzro;
-    if (neginf >= zero) {
+    if (!Risinf(neginf) || neginf >= zero) {
         return_value = 0;
         return return_value;
     }
@@ -73,19 +91,19 @@ INTEGER iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
     }
     //
     posinf = one / newzro;
-    if (posinf <= one) {
+    if (!Risinf(posinf) || posinf <= zero) {
         return_value = 0;
         return return_value;
     }
     //
     neginf = neginf * posinf;
-    if (neginf >= zero) {
+    if (!Risinf(neginf) || neginf >= zero) {
         return_value = 0;
         return return_value;
     }
     //
     posinf = posinf * posinf;
-    if (posinf <= one) {
+    if (!Risinf(posinf) || posinf <= zero) {
         return_value = 0;
         return return_value;
     }
@@ -108,32 +126,32 @@ INTEGER iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
     //
     REAL nan6 = nan5 * zero;
     //
-    if (nan1 == nan1) {
+    if (!Risnan(nan1)) {
         return_value = 0;
         return return_value;
     }
     //
-    if (nan2 == nan2) {
+    if (!Risnan(nan2)) {
         return_value = 0;
         return return_value;
     }
     //
-    if (nan3 == nan3) {
+    if (!Risnan(nan3)) {
         return_value = 0;
         return return_value;
     }
     //
-    if (nan4 == nan4) {
+    if (!Risnan(nan4)) {
         return_value = 0;
         return return_value;
     }
     //
-    if (nan5 == nan5) {
+    if (!Risnan(nan5)) {
         return_value = 0;
         return return_value;
     }
     //
-    if (nan6 == nan6) {
+    if (!Risnan(nan6)) {
         return_value = 0;
         return return_value;
     }
