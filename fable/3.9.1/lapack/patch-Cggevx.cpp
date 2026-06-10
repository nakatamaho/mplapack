--- a/mplapack/reference/Cggevx.cpp
+++ b/mplapack/reference/Cggevx.cpp
@@ -379,7 +379,7 @@
         for (jc = 1; jc <= n; jc = jc + 1) {
             temp = zero;
             for (jr = 1; jr <= n; jr = jr + 1) {
-                temp = max(temp, abs1(vl[(jr - 1) + (jc - 1) * ldvl]));
+                temp = max(temp, cabs1(vl[(jr - 1) + (jc - 1) * ldvl]));
             }
             if (temp < smlnum) {
                 goto statement_50;
@@ -397,7 +397,7 @@
         for (jc = 1; jc <= n; jc = jc + 1) {
             temp = zero;
             for (jr = 1; jr <= n; jr = jr + 1) {
-                temp = max(temp, abs1(vr[(jr - 1) + (jc - 1) * ldvr]));
+                temp = max(temp, cabs1(vr[(jr - 1) + (jc - 1) * ldvr]));
             }
             if (temp < smlnum) {
                 goto statement_80;
