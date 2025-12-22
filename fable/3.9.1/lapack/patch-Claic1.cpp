diff --git a/mplapack/reference/Claic1.cpp b/mplapack/reference/Claic1.cpp
index 3aacf72f..1824a273 100644
--- a/mplapack/reference/Claic1.cpp
+++ b/mplapack/reference/Claic1.cpp
@@ -70,7 +70,7 @@ void Claic1(INTEGER const job, INTEGER const j, COMPLEX *x, REAL const &sest, CO
             } else {
                 s = alpha / s1;
                 c = gamma / s1;
-                tmp = sqrt((s * conj(s) + c * conj(c)).real());
+                tmp = sqrt(s * conj(s) + c * conj(c));
                 s = s / tmp;
                 c = c / tmp;
                 sestpr = s1 * tmp;
@@ -126,12 +126,12 @@ void Claic1(INTEGER const job, INTEGER const j, COMPLEX *x, REAL const &sest, CO
             if (b > zero) {
                 t = (c / (b + sqrt(b * b + c))).real();
             } else {
-                t = (sqrt(b * b + c) - b).real();
+                t = sqrt(b * b + c) - b;
             }
             //
             sine = -(alpha / absest) / t;
             cosine = -(gamma / absest) / (one + t);
-            tmp = sqrt((sine * conj(sine) + cosine * conj(cosine)).real());
+            tmp = sqrt(sine * conj(sine) + cosine * conj(cosine));
             s = sine / tmp;
             c = cosine / tmp;
             sestpr = sqrt(t + one) * absest;
@@ -156,7 +156,7 @@ void Claic1(INTEGER const job, INTEGER const j, COMPLEX *x, REAL const &sest, CO
             s1 = max(abs(sine), abs(cosine));
             s = sine / s1;
             c = cosine / s1;
-            tmp = sqrt((s * conj(s) + c * conj(c)).real());
+            tmp = sqrt(s * conj(s) + c * conj(c));
             s = s / tmp;
             c = c / tmp;
             return;
@@ -224,15 +224,15 @@ void Claic1(INTEGER const job, INTEGER const j, COMPLEX *x, REAL const &sest, CO
                 b = (zeta2 * zeta2 + zeta1 * zeta1 - one) * half;
                 c = zeta1 * zeta1;
                 if (b >= zero) {
-                    t = (-c / (b + sqrt(b * b + c))).real();
+                    t = -c / (b + sqrt(b * b + c));
                 } else {
-                    t = (b - sqrt(b * b + c)).real();
+                    t = b - sqrt(b * b + c);
                 }
                 sine = -(alpha / absest) / t;
                 cosine = -(gamma / absest) / (one + t);
                 sestpr = sqrt(one + t + four * eps * eps * norma) * absest;
             }
-            tmp = sqrt((sine * conj(sine) + cosine * conj(cosine)).real());
+            tmp = sqrt(sine * conj(sine) + cosine * conj(cosine));
             s = sine / tmp;
             c = cosine / tmp;
             return;
