--- a/mplapack/reference/Chgeqz.cpp
+++ b/mplapack/reference/Chgeqz.cpp
@@ -274,7 +274,7 @@
         if (ilast == ilo) {
             goto statement_60;
         } else {
-            if (abs1(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) <= max(safmin, ulp * (abs1(h[(ilast - 1) + (ilast - 1) * ldh]) + abs1(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh])))) {
+            if (cabs1(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) <= max(safmin, ulp * (cabs1(h[(ilast - 1) + (ilast - 1) * ldh]) + cabs1(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh])))) {
                 h[(ilast - 1) + ((ilast - 1) - 1) * ldh] = czero;
                 goto statement_60;
             }
@@ -294,7 +294,7 @@
             if (j == ilo) {
                 ilazro = true;
             } else {
-                if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) <= max(safmin, ulp * (abs1(h[(j - 1) + (j - 1) * ldh]) + abs1(h[((j - 1) - 1) + ((j - 1) - 1) * ldh])))) {
+                if (cabs1(h[(j - 1) + ((j - 1) - 1) * ldh]) <= max(safmin, ulp * (cabs1(h[(j - 1) + (j - 1) * ldh]) + cabs1(h[((j - 1) - 1) + ((j - 1) - 1) * ldh])))) {
                     h[(j - 1) + ((j - 1) - 1) * ldh] = czero;
                     ilazro = true;
                 } else {
@@ -315,7 +315,7 @@
                 //
                 ilazr2 = false;
                 if (!ilazro) {
-                    if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * (ascale * abs1(h[((j + 1) - 1) + (j - 1) * ldh])) <= abs1(h[(j - 1) + (j - 1) * ldh]) * (ascale * atol)) {
+                    if (cabs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * (ascale * cabs1(h[((j + 1) - 1) + (j - 1) * ldh])) <= cabs1(h[(j - 1) + (j - 1) * ldh]) * (ascale * atol)) {
                         ilazr2 = true;
                     }
                 }
@@ -340,7 +340,7 @@
                             h[(jch - 1) + ((jch - 1) - 1) * ldh] = h[(jch - 1) + ((jch - 1) - 1) * ldh] * c;
                         }
                         ilazr2 = false;
-                        if (abs1(t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt]) >= btol) {
+                        if (cabs1(t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt]) >= btol) {
                             if (jch + 1 >= ilast) {
                                 goto statement_60;
                             } else {
@@ -485,11 +485,11 @@
             //
             shift = abi22;
             ctemp = sqrt(abi12) * sqrt(ad21);
-            temp = abs1(ctemp);
+            temp = cabs1(ctemp);
             if (ctemp != zero) {
                 x = half * (ad11 - shift);
-                temp2 = abs1(x);
-                temp = max(temp, abs1(x));
+                temp2 = cabs1(x);
+                temp = max(temp, cabs1(x));
                 y = temp * sqrt(pow2((x / temp)) + pow2((ctemp / temp)));
                 if (temp2 > zero) {
                     if ((x / temp2).real() * y.real() + (x / temp2).imag() * y.imag() < zero) {
@@ -502,7 +502,7 @@
             //
             // Exceptional shift.  Chosen for no particularly good reason.
             //
-            if ((iiter / 20) * 20 == iiter && bscale * abs1(t[(ilast - 1) + (ilast - 1) * ldt]) > safmin) {
+            if ((iiter / 20) * 20 == iiter && bscale * cabs1(t[(ilast - 1) + (ilast - 1) * ldt]) > safmin) {
                 eshift += (ascale * h[(ilast - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
             } else {
                 eshift += (ascale * h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
@@ -515,14 +515,14 @@
         for (j = ilast - 1; j >= ifirst + 1; j = j - 1) {
             istart = j;
             ctemp = ascale * h[(j - 1) + (j - 1) * ldh] - shift * (bscale * t[(j - 1) + (j - 1) * ldt]);
-            temp = abs1(ctemp);
-            temp2 = ascale * abs1(h[((j + 1) - 1) + (j - 1) * ldh]);
+            temp = cabs1(ctemp);
+            temp2 = ascale * cabs1(h[((j + 1) - 1) + (j - 1) * ldh]);
             tempr = max(temp, temp2);
             if (tempr < one && tempr != zero) {
                 temp = temp / tempr;
                 temp2 = temp2 / tempr;
             }
-            if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * temp2 <= temp * atol) {
+            if (cabs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * temp2 <= temp * atol) {
                 goto statement_90;
             }
         }
