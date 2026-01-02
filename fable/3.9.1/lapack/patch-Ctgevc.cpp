--- a/mplapack/reference/Ctgevc.cpp
+++ b/mplapack/reference/Ctgevc.cpp
@@ -195,19 +195,19 @@
     // part of A and B to check for possible overflow in the triangular
     // solver.
     //
-    anorm = abs1(s[0]);
-    bnorm = abs1(p[0]);
+    anorm = cabs1(s[0]);
+    bnorm = cabs1(p[0]);
     rwork[1 - 1] = zero;
     rwork[(n + 1) - 1] = zero;
     for (j = 2; j <= n; j = j + 1) {
         rwork[j - 1] = zero;
         rwork[(n + j) - 1] = zero;
         for (i = 1; i <= j - 1; i = i + 1) {
-            rwork[j - 1] += abs1(s[(i - 1) + (j - 1) * lds]);
-            rwork[(n + j) - 1] += abs1(p[(i - 1) + (j - 1) * ldp]);
+            rwork[j - 1] += cabs1(s[(i - 1) + (j - 1) * lds]);
+            rwork[(n + j) - 1] += cabs1(p[(i - 1) + (j - 1) * ldp]);
         }
-        anorm = max(anorm, rwork[j - 1] + abs1(s[(j - 1) + (j - 1) * lds]));
-        bnorm = max(bnorm, rwork[(n + j) - 1] + abs1(p[(j - 1) + (j - 1) * ldp]));
+        anorm = max(anorm, rwork[j - 1] + cabs1(s[(j - 1) + (j - 1) * lds]));
+        bnorm = max(bnorm, rwork[(n + j) - 1] + cabs1(p[(j - 1) + (j - 1) * ldp]));
     }
     //
     ascale = one / max(anorm, safmin);
@@ -229,7 +229,7 @@
             if (ilcomp) {
                 ieig++;
                 //
-                if (abs1(s[(je - 1) + (je - 1) * lds]) <= safmin && abs(p[(je - 1) + (je - 1) * ldp].real()) <= safmin) {
+                if (cabs1(s[(je - 1) + (je - 1) * lds]) <= safmin && abs(p[(je - 1) + (je - 1) * ldp].real()) <= safmin) {
                     //
                     // Singular matrix pencil -- return unit eigenvector
                     //
@@ -245,7 +245,7 @@
                 // H
                 // y  ( a A - b B ) = 0
                 //
-                temp = one / max(abs1(s[(je - 1) + (je - 1) * lds]) * ascale, abs(p[(je - 1) + (je - 1) * ldp].real()) * bscale, safmin);
+                temp = one / max(cabs1(s[(je - 1) + (je - 1) * lds]) * ascale, abs(p[(je - 1) + (je - 1) * ldp].real()) * bscale, safmin);
                 salpha = (temp * s[(je - 1) + (je - 1) * lds]) * ascale;
                 sbeta = (temp * p[(je - 1) + (je - 1) * ldp].real()) * bscale;
                 acoeff = sbeta * ascale;
@@ -254,17 +254,17 @@
                 // Scale to avoid underflow
                 //
                 lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
-                lsb = abs1(salpha) >= safmin && abs1(bcoeff) < small;
+                lsb = cabs1(salpha) >= safmin && cabs1(bcoeff) < small;
                 //
                 scale = one;
                 if (lsa) {
                     scale = (small / abs(sbeta)) * min(anorm, big);
                 }
                 if (lsb) {
-                    scale = max(scale, (small / abs1(salpha)) * min(bnorm, big));
+                    scale = max(scale, (small / cabs1(salpha)) * min(bnorm, big));
                 }
                 if (lsa || lsb) {
-                    scale = min(scale, one / (safmin * max(one, abs(acoeff), abs1(bcoeff))));
+                    scale = min(scale, one / (safmin * max(one, abs(acoeff), cabs1(bcoeff))));
                     if (lsa) {
                         acoeff = ascale * (scale * sbeta);
                     } else {
@@ -278,7 +278,7 @@
                 }
                 //
                 acoefa = abs(acoeff);
-                bcoefa = abs1(bcoeff);
+                bcoefa = cabs1(bcoeff);
                 xmax = one;
                 for (jr = 1; jr <= n; jr = jr + 1) {
                     work[jr - 1] = czero;
@@ -321,13 +321,13 @@
                     // with scaling and perturbation of the denominator
                     //
                     d = conj(acoeff * s[(j - 1) + (j - 1) * lds] - bcoeff * p[(j - 1) + (j - 1) * ldp]);
-                    if (abs1(d) <= dmin) {
+                    if (cabs1(d) <= dmin) {
                         d = COMPLEX(dmin);
                     }
                     //
-                    if (abs1(d) < one) {
-                        if (abs1(sum) >= bignum * abs1(d)) {
-                            temp = one / abs1(sum);
+                    if (cabs1(d) < one) {
+                        if (cabs1(sum) >= bignum * cabs1(d)) {
+                            temp = one / cabs1(sum);
                             for (jr = je; jr <= j - 1; jr = jr + 1) {
                                 work[jr - 1] = temp * work[jr - 1];
                             }
@@ -336,7 +336,7 @@
                         }
                     }
                     work[j - 1] = Cladiv(-sum, d);
-                    xmax = max(xmax, abs1(work[j - 1]));
+                    xmax = max(xmax, cabs1(work[j - 1]));
                 }
                 //
                 // Back transform eigenvector if HOWMNY='B'.
@@ -354,7 +354,7 @@
                 //
                 xmax = zero;
                 for (jr = ibeg; jr <= n; jr = jr + 1) {
-                    xmax = max(xmax, abs1(work[((isrc - 1) * n + jr) - 1]));
+                    xmax = max(xmax, cabs1(work[((isrc - 1) * n + jr) - 1]));
                 }
                 //
                 if (xmax > safmin) {
@@ -391,7 +391,7 @@
             if (ilcomp) {
                 ieig = ieig - 1;
                 //
-                if (abs1(s[(je - 1) + (je - 1) * lds]) <= safmin && abs(p[(je - 1) + (je - 1) * ldp].real()) <= safmin) {
+                if (cabs1(s[(je - 1) + (je - 1) * lds]) <= safmin && abs(p[(je - 1) + (je - 1) * ldp].real()) <= safmin) {
                     //
                     // Singular matrix pencil -- return unit eigenvector
                     //
@@ -407,7 +407,7 @@
                 //
                 // ( a A - b B ) x  = 0
                 //
-                temp = one / max(abs1(s[(je - 1) + (je - 1) * lds]) * ascale, abs(p[(je - 1) + (je - 1) * ldp].real()) * bscale, safmin);
+                temp = one / max(cabs1(s[(je - 1) + (je - 1) * lds]) * ascale, abs(p[(je - 1) + (je - 1) * ldp].real()) * bscale, safmin);
                 salpha = (temp * s[(je - 1) + (je - 1) * lds]) * ascale;
                 sbeta = (temp * p[(je - 1) + (je - 1) * ldp].real()) * bscale;
                 acoeff = sbeta * ascale;
@@ -416,17 +416,17 @@
                 // Scale to avoid underflow
                 //
                 lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
-                lsb = abs1(salpha) >= safmin && abs1(bcoeff) < small;
+                lsb = cabs1(salpha) >= safmin && cabs1(bcoeff) < small;
                 //
                 scale = one;
                 if (lsa) {
                     scale = (small / abs(sbeta)) * min(anorm, big);
                 }
                 if (lsb) {
-                    scale = max(scale, (small / abs1(salpha)) * min(bnorm, big));
+                    scale = max(scale, (small / cabs1(salpha)) * min(bnorm, big));
                 }
                 if (lsa || lsb) {
-                    scale = min(scale, one / (safmin * max(one, abs(acoeff), abs1(bcoeff))));
+                    scale = min(scale, one / (safmin * max(one, abs(acoeff), cabs1(bcoeff))));
                     if (lsa) {
                         acoeff = ascale * (scale * sbeta);
                     } else {
@@ -440,7 +440,7 @@
                 }
                 //
                 acoefa = abs(acoeff);
-                bcoefa = abs1(bcoeff);
+                bcoefa = cabs1(bcoeff);
                 xmax = one;
                 for (jr = 1; jr <= n; jr = jr + 1) {
                     work[jr - 1] = czero;
@@ -464,13 +464,13 @@
                     // with scaling and perturbation of the denominator
                     //
                     d = acoeff * s[(j - 1) + (j - 1) * lds] - bcoeff * p[(j - 1) + (j - 1) * ldp];
-                    if (abs1(d) <= dmin) {
+                    if (cabs1(d) <= dmin) {
                         d = COMPLEX(dmin);
                     }
                     //
-                    if (abs1(d) < one) {
-                        if (abs1(work[j - 1]) >= bignum * abs1(d)) {
-                            temp = one / abs1(work[j - 1]);
+                    if (cabs1(d) < one) {
+                        if (cabs1(work[j - 1]) >= bignum * cabs1(d)) {
+                            temp = one / cabs1(work[j - 1]);
                             for (jr = 1; jr <= je; jr = jr + 1) {
                                 work[jr - 1] = temp * work[jr - 1];
                             }
@@ -483,8 +483,8 @@
                         //
                         // w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
                         //
-                        if (abs1(work[j - 1]) > one) {
-                            temp = one / abs1(work[j - 1]);
+                        if (cabs1(work[j - 1]) > one) {
+                            temp = one / cabs1(work[j - 1]);
                             if (acoefa * rwork[j - 1] + bcoefa * rwork[(n + j) - 1] >= bignum * temp) {
                                 for (jr = 1; jr <= je; jr = jr + 1) {
                                     work[jr - 1] = temp * work[jr - 1];
@@ -515,7 +515,7 @@
                 //
                 xmax = zero;
                 for (jr = 1; jr <= iend; jr = jr + 1) {
-                    xmax = max(xmax, abs1(work[((isrc - 1) * n + jr) - 1]));
+                    xmax = max(xmax, cabs1(work[((isrc - 1) * n + jr) - 1]));
                 }
                 //
                 if (xmax > safmin) {
