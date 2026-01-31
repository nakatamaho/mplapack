--- Rckcsd.cpp
+++ Rckcsd.cpp
@@ -108,9 +108,11 @@ void Rckcsd(INTEGER const nm, INTEGER *mval, INTEGER *pval, INTEGER *qval, INTEG
     INTEGER iinfo = 0;
     INTEGER r = 0;
     INTEGER i = 0;
-    const REAL piover2 = 1.5707963267948966192313216916397514421;
+    const REAL dummy = 0.0;
+    const REAL two = 2.0;
+    const REAL piover2 = pi(dummy) / 2.0;
     INTEGER j = 0;
-    const REAL orth = 0.000000000001;
+    const REAL orth = 1.0e-12;
     const REAL ten = 10.0;
     const REAL gapdigit = 18.0;
     const REAL zero = 0.0;
