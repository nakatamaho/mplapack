diff --git a/mplapack/reference/Clartg.cpp b/mplapack/reference/Clartg.cpp
index c42acaff..53c7325f 100644
--- a/mplapack/reference/Clartg.cpp
+++ b/mplapack/reference/Clartg.cpp
@@ -51,12 +51,13 @@ void Clartg(COMPLEX const &f, COMPLEX const &g, REAL &cs, COMPLEX &sn, COMPLEX &
     REAL dr = 0.0;
     REAL di = 0.0;
     INTEGER i = 0;
+    abs1(ff) = max(abs(ff.real()), abs(ff.imag()));
     //
     safmin = Rlamch("S");
     eps = Rlamch("E");
     safmn2 = pow(Rlamch("B"), castINTEGER(log(safmin / eps) / log(Rlamch("B")) / two));
     safmx2 = one / safmn2;
-    scale = max(cabs_inf(f), cabs_inf(g));
+    scale = max(abs1(f), abs1(g));
     fs = f;
     gs = g;
     count = 0;
@@ -113,7 +114,7 @@ void Clartg(COMPLEX const &f, COMPLEX const &g, REAL &cs, COMPLEX &sn, COMPLEX &
         cs = f2s / g2s;
         // Make sure abs(FF) = 1
         // Do complex/real division explicitly with 2 real divisions
-        if (cabs_inf(f) > one) {
+        if (abs1(f) > one) {
             d = Rlapy2(f.real(), f.imag());
             ff = COMPLEX(f.real() / d, f.imag() / d);
         } else {
