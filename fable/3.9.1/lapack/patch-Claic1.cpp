--- a/mplapack/reference/Claic1.cpp
+++ b/mplapack/reference/Claic1.cpp
@@ -77,7 +77,7 @@
             } else {
                 s = alpha / s1;
                 c = gamma / s1;
-                tmp = sqrt(s * conj(s) + c * conj(c));
+                tmp = sqrt((s * conj(s) + c * conj(c)).real());
                 s = s / tmp;
                 c = c / tmp;
                 sestpr = s1 * tmp;
@@ -133,12 +133,12 @@
             if (b > zero) {
                 t = (c / (b + sqrt(b * b + c))).real();
             } else {
-                t = sqrt(b * b + c) - b;
+                t = (sqrt(b * b + c) - b).real();
             }
             //
             sine = -(alpha / absest) / t;
             cosine = -(gamma / absest) / (one + t);
-            tmp = sqrt(sine * conj(sine) + cosine * conj(cosine));
+            tmp = sqrt((sine * conj(sine) + cosine * conj(cosine)).real());
             s = sine / tmp;
             c = cosine / tmp;
             sestpr = sqrt(t + one) * absest;
@@ -163,7 +163,7 @@
             s1 = max(abs(sine), abs(cosine));
             s = sine / s1;
             c = cosine / s1;
-            tmp = sqrt(s * conj(s) + c * conj(c));
+            tmp = sqrt((s * conj(s) + c * conj(c)).real());
             s = s / tmp;
             c = c / tmp;
             return;
@@ -231,15 +231,15 @@
                 b = (zeta2 * zeta2 + zeta1 * zeta1 - one) * half;
                 c = zeta1 * zeta1;
                 if (b >= zero) {
-                    t = -c / (b + sqrt(b * b + c));
+                    t = (-c / (b + sqrt(b * b + c))).real();
                 } else {
-                    t = b - sqrt(b * b + c);
+                    t = (b - sqrt(b * b + c)).real();
                 }
                 sine = -(alpha / absest) / t;
                 cosine = -(gamma / absest) / (one + t);
                 sestpr = sqrt(one + t + four * eps * eps * norma) * absest;
             }
-            tmp = sqrt(sine * conj(sine) + cosine * conj(cosine));
+            tmp = sqrt((sine * conj(sine) + cosine * conj(cosine)).real());
             s = sine / tmp;
             c = cosine / tmp;
             return;
