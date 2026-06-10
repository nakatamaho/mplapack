--- a/mplapack/reference/Clartg.cpp
+++ b/mplapack/reference/Clartg.cpp
@@ -36,6 +36,13 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL cabs_inf(COMPLEX ff) { return max(abs(ff.real()), abs(ff.imag())); }
+inline REAL cabssq(COMPLEX ff) {
+    REAL temp;
+    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
+    return temp;
+}
+
 void Clartg(COMPLEX const f, COMPLEX const g, REAL &cs, COMPLEX &sn, COMPLEX &r) {
     COMPLEX ff = 0.0;
     REAL safmin = 0.0;
@@ -58,13 +65,12 @@
     REAL dr = 0.0;
     REAL di = 0.0;
     INTEGER i = 0;
-    abs1(ff) = max(abs(ff.real()), abs(ff.imag()));
     //
     safmin = Rlamch("S");
     eps = Rlamch("E");
     safmn2 = pow(Rlamch("B"), castINTEGER(log(safmin / eps) / log(Rlamch("B")) / two));
     safmx2 = one / safmn2;
-    scale = max(abs1(f), abs1(g));
+    scale = max(cabs_inf(f), cabs_inf(g));
     fs = f;
     gs = g;
     count = 0;
@@ -93,8 +99,8 @@
             goto statement_20;
         }
     }
-    f2 = abssq(fs);
-    g2 = abssq(gs);
+    f2 = cabssq(fs);
+    g2 = cabssq(gs);
     if (f2 <= max(g2, one) * safmin) {
         //
         // This is a rare case: F is very small.
@@ -121,7 +127,7 @@
         cs = f2s / g2s;
         // Make sure abs(FF) = 1
         // Do complex/real division explicitly with 2 real divisions
-        if (abs1(f) > one) {
+        if (cabs_inf(f) > one) {
             d = Rlapy2(f.real(), f.imag());
             ff = COMPLEX(f.real() / d, f.imag() / d);
         } else {
