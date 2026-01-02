--- a/mplapack/reference/Clags2.cpp
+++ b/mplapack/reference/Clags2.cpp
@@ -113,16 +113,16 @@
             vb11r = csr * b1;
             vb12 = csr * b2 + d1 * snr * b3;
             //
-            aua12 = abs(csl) * abs1(a2) + abs(snl) * abs(a3);
-            avb12 = abs(csr) * abs1(b2) + abs(snr) * abs(b3);
+            aua12 = abs(csl) * cabs1(a2) + abs(snl) * abs(a3);
+            avb12 = abs(csr) * cabs1(b2) + abs(snr) * abs(b3);
             //
             // zero (1,2) elements of U**H *A and V**H *B
             //
-            if ((abs(ua11r) + abs1(ua12)) == zero) {
+            if ((abs(ua11r) + cabs1(ua12)) == zero) {
                 Clartg(-COMPLEX(vb11r), conj(vb12), csq, snq, r);
-            } else if ((abs(vb11r) + abs1(vb12)) == zero) {
+            } else if ((abs(vb11r) + cabs1(vb12)) == zero) {
                 Clartg(-COMPLEX(ua11r), conj(ua12), csq, snq, r);
-            } else if (aua12 / (abs(ua11r) + abs1(ua12)) <= avb12 / (abs(vb11r) + abs1(vb12))) {
+            } else if (aua12 / (abs(ua11r) + cabs1(ua12)) <= avb12 / (abs(vb11r) + cabs1(vb12))) {
                 Clartg(-COMPLEX(ua11r), conj(ua12), csq, snq, r);
             } else {
                 Clartg(-COMPLEX(vb11r), conj(vb12), csq, snq, r);
@@ -144,16 +144,16 @@
             vb21 = -conj(d1) * snr * b1;
             vb22 = -conj(d1) * snr * b2 + csr * b3;
             //
-            aua22 = abs(snl) * abs1(a2) + abs(csl) * abs(a3);
-            avb22 = abs(snr) * abs1(b2) + abs(csr) * abs(b3);
+            aua22 = abs(snl) * cabs1(a2) + abs(csl) * abs(a3);
+            avb22 = abs(snr) * cabs1(b2) + abs(csr) * abs(b3);
             //
             // zero (2,2) elements of U**H *A and V**H *B, and then swap.
             //
-            if ((abs1(ua21) + abs1(ua22)) == zero) {
+            if ((cabs1(ua21) + cabs1(ua22)) == zero) {
                 Clartg(-conj(vb21), conj(vb22), csq, snq, r);
-            } else if ((abs1(vb21) + abs(vb22)) == zero) {
+            } else if ((cabs1(vb21) + abs(vb22)) == zero) {
                 Clartg(-conj(ua21), conj(ua22), csq, snq, r);
-            } else if (aua22 / (abs1(ua21) + abs1(ua22)) <= avb22 / (abs1(vb21) + abs1(vb22))) {
+            } else if (aua22 / (cabs1(ua21) + cabs1(ua22)) <= avb22 / (cabs1(vb21) + cabs1(vb22))) {
                 Clartg(-conj(ua21), conj(ua22), csq, snq, r);
             } else {
                 Clartg(-conj(vb21), conj(vb22), csq, snq, r);
@@ -204,16 +204,16 @@
             vb21 = -d1 * snl * b1 + csl * b2;
             vb22r = csl * b3;
             //
-            aua21 = abs(snr) * abs(a1) + abs(csr) * abs1(a2);
-            avb21 = abs(snl) * abs(b1) + abs(csl) * abs1(b2);
+            aua21 = abs(snr) * abs(a1) + abs(csr) * cabs1(a2);
+            avb21 = abs(snl) * abs(b1) + abs(csl) * cabs1(b2);
             //
             // zero (2,1) elements of U**H *A and V**H *B.
             //
-            if ((abs1(ua21) + abs(ua22r)) == zero) {
+            if ((cabs1(ua21) + abs(ua22r)) == zero) {
                 Clartg(COMPLEX(vb22r), vb21, csq, snq, r);
-            } else if ((abs1(vb21) + abs(vb22r)) == zero) {
+            } else if ((cabs1(vb21) + abs(vb22r)) == zero) {
                 Clartg(COMPLEX(ua22r), ua21, csq, snq, r);
-            } else if (aua21 / (abs1(ua21) + abs(ua22r)) <= avb21 / (abs1(vb21) + abs(vb22r))) {
+            } else if (aua21 / (cabs1(ua21) + abs(ua22r)) <= avb21 / (cabs1(vb21) + abs(vb22r))) {
                 Clartg(COMPLEX(ua22r), ua21, csq, snq, r);
             } else {
                 Clartg(COMPLEX(vb22r), vb21, csq, snq, r);
@@ -235,16 +235,16 @@
             vb11 = csl * b1 + conj(d1) * snl * b2;
             vb12 = conj(d1) * snl * b3;
             //
-            aua11 = abs(csr) * abs(a1) + abs(snr) * abs1(a2);
-            avb11 = abs(csl) * abs(b1) + abs(snl) * abs1(b2);
+            aua11 = abs(csr) * abs(a1) + abs(snr) * cabs1(a2);
+            avb11 = abs(csl) * abs(b1) + abs(snl) * cabs1(b2);
             //
             // zero (1,1) elements of U**H *A and V**H *B, and then swap.
             //
-            if ((abs1(ua11) + abs1(ua12)) == zero) {
+            if ((cabs1(ua11) + cabs1(ua12)) == zero) {
                 Clartg(vb12, vb11, csq, snq, r);
-            } else if ((abs1(vb11) + abs1(vb12)) == zero) {
+            } else if ((cabs1(vb11) + cabs1(vb12)) == zero) {
                 Clartg(ua12, ua11, csq, snq, r);
-            } else if (aua11 / (abs1(ua11) + abs1(ua12)) <= avb11 / (abs1(vb11) + abs1(vb12))) {
+            } else if (aua11 / (cabs1(ua11) + cabs1(ua12)) <= avb11 / (cabs1(vb11) + cabs1(vb12))) {
                 Clartg(ua12, ua11, csq, snq, r);
             } else {
                 Clartg(vb12, vb11, csq, snq, r);
