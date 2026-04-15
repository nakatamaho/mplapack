--- a/mplapack/reference/Clartg.cpp	2026-03-21 11:38:06.661888232 +0900
+++ b/mplapack/reference/Clartg.cpp	2026-03-21 11:42:48.544846415 +0900
@@ -35,17 +35,19 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <mplapack_arithmetic_params.h>
 
 void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
     COMPLEX t = 0.0;
-    const REAL safmin = 1.0 + unhandled;
+    const auto &ap = mplapack::get_arithmetic_params<REAL>();
+    const REAL safmin = ap.safmin;
     REAL rtmin = sqrt(safmin);
     //
     const COMPLEX czero = COMPLEX(0.0, 0.0);
     const REAL one = 1.0;
     const REAL zero = 0.0;
     REAL g1 = 0.0;
-    const REAL safmax = 1.0 + unhandled;
+    const REAL safmax = ap.safmax;
     REAL rtmax = 0.0;
     REAL g2 = 0.0;
     REAL d = 0.0;

diff --git a/mplapack/reference/Clartg.cpp b/mplapack/reference/Clartg.cpp
index bc25c968..3d2d8bad 100644
--- a/mplapack/reference/Clartg.cpp
+++ b/mplapack/reference/Clartg.cpp
@@ -37,6 +37,12 @@
 #include <mplapack.h>
 #include <mplapack_arithmetic_params.h>
 
+inline REAL abssq(COMPLEX ff) {
+    REAL temp;
+    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
+    return temp;
+}
+
 void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
     COMPLEX t = 0.0;
     const auto &ap = mplapack::get_arithmetic_params<REAL>();

diff --git a/mplapack/reference/Clartg.cpp b/mplapack/reference/Clartg.cpp
index 3d2d8bad..ce94f6c 100644
--- a/mplapack/reference/Clartg.cpp
+++ b/mplapack/reference/Clartg.cpp
@@ -108,7 +108,16 @@ void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
         f1 = max(abs(f.real()), abs(f.imag()));
         g1 = max(abs(g.real()), abs(g.imag()));
         rtmax = sqrt(safmax / 4);
-        if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
+        REAL fsmall = min(abs(f.real()), abs(f.imag()));
+        REAL gsmall = min(abs(g.real()), abs(g.imag()));
+        REAL split_tol = sqrt(ap.eps);
+        bool use_unscaled = f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax;
+        if (use_unscaled && (f1 < one || g1 < one)) {
+            bool f_has_split = fsmall != zero && fsmall / f1 > split_tol;
+            bool g_has_split = gsmall != zero && gsmall / g1 > split_tol;
+            use_unscaled = !(f_has_split || g_has_split);
+        }
+        if (use_unscaled) {
             //
             // Use unscaled algorithm
             //
@@ -123,7 +132,7 @@ void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
                 rtmax = rtmax * 2;
                 if (f2 > rtmin && h2 < rtmax) {
                     // safmin <= sqrt( f2*h2 ) <= safmax
-                    s = conj(g) * (f / sqrt(f2 * h2));
+                    s = (f / sqrt(f2)) * (conj(g) / sqrt(h2));
                 } else {
                     s = conj(g) * (r / h2);
                 }
@@ -179,7 +188,7 @@ void Clartg(COMPLEX const f, COMPLEX const g, REAL &c, COMPLEX &s, COMPLEX &r) {
                 rtmax = rtmax * 2;
                 if (f2 > rtmin && h2 < rtmax) {
                     // safmin <= sqrt( f2*h2 ) <= safmax
-                    s = conj(gs) * (fs / sqrt(f2 * h2));
+                    s = (fs / sqrt(f2)) * (conj(gs) / sqrt(h2));
                 } else {
                     s = conj(gs) * (r / h2);
                 }
