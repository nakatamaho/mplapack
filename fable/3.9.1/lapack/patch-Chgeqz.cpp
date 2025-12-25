diff --git a/mplapack/reference/Chgeqz.cpp b/mplapack/reference/Chgeqz.cpp
index 98595357..452d156b 100644
--- a/mplapack/reference/Chgeqz.cpp
+++ b/mplapack/reference/Chgeqz.cpp
@@ -274,7 +274,7 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
         if (ilast == ilo) {
             goto statement_60;
         } else {
-            if (cabs1(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) <= max(safmin, ulp * (cabs1(h[(ilast - 1) + (ilast - 1) * ldh]) + cabs1(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh])))) {
+            if (abs1(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) <= max(safmin, ulp * (abs1(h[(ilast - 1) + (ilast - 1) * ldh]) + abs1(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh])))) {
                 h[(ilast - 1) + ((ilast - 1) - 1) * ldh] = czero;
                 goto statement_60;
             }
@@ -294,7 +294,7 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
             if (j == ilo) {
                 ilazro = true;
             } else {
-                if (cabs1(h[(j - 1) + ((j - 1) - 1) * ldh]) <= max(safmin, ulp * (cabs1(h[(j - 1) + (j - 1) * ldh]) + cabs1(h[((j - 1) - 1) + ((j - 1) - 1) * ldh])))) {
+                if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) <= max(safmin, ulp * (abs1(h[(j - 1) + (j - 1) * ldh]) + abs1(h[((j - 1) - 1) + ((j - 1) - 1) * ldh])))) {
                     h[(j - 1) + ((j - 1) - 1) * ldh] = czero;
                     ilazro = true;
                 } else {
@@ -315,7 +315,7 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
                 //
                 ilazr2 = false;
                 if (!ilazro) {
-                    if (cabs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * (ascale * cabs1(h[((j + 1) - 1) + (j - 1) * ldh])) <= cabs1(h[(j - 1) + (j - 1) * ldh]) * (ascale * atol)) {
+                    if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * (ascale * abs1(h[((j + 1) - 1) + (j - 1) * ldh])) <= abs1(h[(j - 1) + (j - 1) * ldh]) * (ascale * atol)) {
                         ilazr2 = true;
                     }
                 }
@@ -340,7 +340,7 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
                             h[(jch - 1) + ((jch - 1) - 1) * ldh] = h[(jch - 1) + ((jch - 1) - 1) * ldh] * c;
                         }
                         ilazr2 = false;
-                        if (cabs1(t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt]) >= btol) {
+                        if (abs1(t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt]) >= btol) {
                             if (jch + 1 >= ilast) {
                                 goto statement_60;
                             } else {
@@ -485,11 +485,11 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
             //
             shift = abi22;
             ctemp = sqrt(abi12) * sqrt(ad21);
-            temp = cabs1(ctemp);
+            temp = abs1(ctemp);
             if (ctemp != zero) {
                 x = half * (ad11 - shift);
-                temp2 = cabs1(x);
-                temp = max(temp, cabs1(x));
+                temp2 = abs1(x);
+                temp = max(temp, abs1(x));
                 y = temp * sqrt(pow2((x / temp)) + pow2((ctemp / temp)));
                 if (temp2 > zero) {
                     if ((x / temp2).real() * y.real() + (x / temp2).imag() * y.imag() < zero) {
@@ -502,7 +502,7 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
             //
             // Exceptional shift.  Chosen for no particularly good reason.
             //
-            if ((iiter / 20) * 20 == iiter && bscale * cabs1(t[(ilast - 1) + (ilast - 1) * ldt]) > safmin) {
+            if ((iiter / 20) * 20 == iiter && bscale * abs1(t[(ilast - 1) + (ilast - 1) * ldt]) > safmin) {
                 eshift += (ascale * h[(ilast - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
             } else {
                 eshift += (ascale * h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
@@ -515,14 +515,14 @@ void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const
         for (j = ilast - 1; j >= ifirst + 1; j = j - 1) {
             istart = j;
             ctemp = ascale * h[(j - 1) + (j - 1) * ldh] - shift * (bscale * t[(j - 1) + (j - 1) * ldt]);
-            temp = cabs1(ctemp);
-            temp2 = ascale * cabs1(h[((j + 1) - 1) + (j - 1) * ldh]);
+            temp = abs1(ctemp);
+            temp2 = ascale * abs1(h[((j + 1) - 1) + (j - 1) * ldh]);
             tempr = max(temp, temp2);
             if (tempr < one && tempr != zero) {
                 temp = temp / tempr;
                 temp2 = temp2 / tempr;
             }
-            if (cabs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * temp2 <= temp * atol) {
+            if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * temp2 <= temp * atol) {
                 goto statement_90;
             }
         }
diff --git a/mplapack/reference/Clahef_aa.cpp b/mplapack/reference/Clahef_aa.cpp
index 58423747..14540bcc 100644
--- a/mplapack/reference/Clahef_aa.cpp
+++ b/mplapack/reference/Clahef_aa.cpp
@@ -131,7 +131,7 @@ void Clahef_aa(const char *uplo, INTEGER const j1, INTEGER const m, INTEGER cons
             //
             // Apply hermitian pivot
             //
-            if ((i2 != 2) && (piv != zero)) {
+            if ((i2 != 2) && (piv != 0)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
@@ -279,7 +279,7 @@ void Clahef_aa(const char *uplo, INTEGER const j1, INTEGER const m, INTEGER cons
             //
             // Apply hermitian pivot
             //
-            if ((i2 != 2) && (piv != zero)) {
+            if ((i2 != 2) && (piv != 0)) {
                 //
                 // Swap WORK(I1) and WORK(I2)
                 //
