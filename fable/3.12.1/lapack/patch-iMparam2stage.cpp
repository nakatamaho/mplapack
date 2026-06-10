--- a/mplapack/reference/iMparam2stage.cpp
+++ b/mplapack/reference/iMparam2stage.cpp
@@ -36,8 +36,12 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-INTEGER
-iMparam2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const ni, INTEGER const nbi, INTEGER const ibi, INTEGER const nxi) {
+#include <ctype.h>
+#include <string.h>
+
+#define SUBNAM_LEN 128
+
+INTEGER iMparam2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const ni, INTEGER const nbi, INTEGER const ibi, INTEGER const nxi) {
     INTEGER return_value = 0;
     //
     // Invalid value for ISPEC
@@ -50,79 +54,43 @@
     // Get the number of threads
     //
     INTEGER nthreads = 1;
-    // WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC
     //
-    char subnam[12];
-    INTEGER ic = 0;
-    INTEGER iz = 0;
+    char subnam[SUBNAM_LEN];
     INTEGER i = 0;
+    INTEGER name_len = 16;
     char prec;
-    char algo[3];
-    char stag[5];
+    char algo[4];
+    char stag[6];
     bool rprec = false;
     bool cprec = false;
+
+    memset(subnam, 0, SUBNAM_LEN);
+    strncpy(subnam, name, name_len);
+    subnam[name_len] = '\0';
+    for (i = 0; i < name_len; i++) {
+        subnam[i] = toupper((unsigned char)subnam[i]);
+    }
+
     if (ispec != 19) {
         //
         // Convert NAME to upper case if the first character is lower case.
         //
         return_value = -1;
-        subnam = *name;
-        ic = fem::ichar(subnam[0]);
-        iz = fem::ichar("Z");
-        if (iz == 90 || iz == 122) {
-            //
-            // ASCII character set
-            //
-            if (ic >= 97 && ic <= 122) {
-                subnam[0] = fem::fchar(ic - 32);
-                for (i = 2; i <= 12; i = i + 1) {
-                    ic = fem::ichar(subnam(i, i));
-                    if (ic >= 97 && ic <= 122) {
-                        subnam(i, i) = fem::fchar(ic - 32);
-                    }
-                }
-            }
-            //
-        } else if (iz == 233 || iz == 169) {
-            //
-            // EBCDIC character set
-            //
-            if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
-                subnam[0] = fem::fchar(ic + 64);
-                for (i = 2; i <= 12; i = i + 1) {
-                    ic = fem::ichar(subnam(i, i));
-                    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
-                        subnam(i, i) = fem::fchar(ic + 64);
-                    }
-                }
-            }
-            //
-        } else if (iz == 218 || iz == 250) {
-            //
-            // Prime machines:  ASCII+128
-            //
-            if (ic >= 225 && ic <= 250) {
-                subnam[0] = fem::fchar(ic - 32);
-                for (i = 2; i <= 12; i = i + 1) {
-                    ic = fem::ichar(subnam(i, i));
-                    if (ic >= 225 && ic <= 250) {
-                        subnam(i, i) = fem::fchar(ic - 32);
-                    }
-                }
-            }
-        }
         //
         prec = subnam[0];
-        algo[0] = subnam[(4 - 1)];
-        algo[1] = subnam[(5 - 1)];
-        algo[2] = subnam[(6 - 1)];
-        stag[0] = subnam[(8 - 1)];
-        stag[1] = subnam[(9 - 1)];
-        stag[2] = subnam[(10 - 1)];
-        stag[3] = subnam[(11 - 1)];
-        stag[4] = subnam[(12 - 1)];
-        rprec = prec == "S" || prec == "D";
-        cprec = prec == "C" || prec == "Z";
+        algo[0] = subnam[3];
+        algo[1] = subnam[4];
+        algo[2] = subnam[5];
+        algo[3] = '\0';
+
+        stag[0] = subnam[7];
+        stag[1] = subnam[8];
+        stag[2] = subnam[9];
+        stag[3] = subnam[10];
+        stag[4] = subnam[11];
+        stag[5] = '\0';
+        rprec = (prec == 'R');
+        cprec = (prec == 'C');
         //
         // Invalid value for PRECISION
         //
@@ -131,8 +99,6 @@
             return return_value;
         }
     }
-    // WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC,
-    // $           '   ALGO ',ALGO,'    STAGE ',STAG
     //
     INTEGER kd = 0;
     INTEGER ib = 0;
@@ -217,26 +183,36 @@
         // + (KD+1)*N
         lwork = -1;
         subnam[0] = prec;
-        subnam(2, 6) = "GEQRF";
+        subnam[1] = 'G';
+        subnam[2] = 'E';
+        subnam[3] = 'Q';
+        subnam[4] = 'R';
+        subnam[5] = 'F';
+        subnam[6] = '\0';
         qroptnb = iMlaenv(1, subnam, " ", ni, nbi, -1, -1);
-        subnam(2, 6) = "GELQF";
+        subnam[1] = 'G';
+        subnam[2] = 'E';
+        subnam[3] = 'L';
+        subnam[4] = 'Q';
+        subnam[5] = 'F';
+        subnam[6] = '\0';
         lqoptnb = iMlaenv(1, subnam, " ", nbi, ni, -1, -1);
         // Could be QR or LQ for TRD and the max for BRD
         factoptnb = max(qroptnb, lqoptnb);
-        if (algo == "TRD") {
-            if (stag == "2STAG") {
-                lwork = ni * nbi + ni * max(nbi + 1, factoptnb) + max(2 * nbi * nbi, nbi * nthreads) + (nbi + 1) * ni;
-            } else if ((stag == "HE2HB") || (stag == "SY2SB")) {
+        if (strncmp(algo, "TRD", 3) == 0) {
+            if (strncmp(stag, "2STAG", 5) == 0) {
+                lwork = ni * nbi + ni * max(nbi + 1, factoptnb) + max((INTEGER)2 * nbi * nbi, nbi * nthreads) + (nbi + 1) * ni;
+            } else if ((strncmp(stag, "HE2HB", 5) == 0) || (strncmp(stag, "SY2SB", 5) == 0)) {
                 lwork = ni * nbi + ni * max(nbi, factoptnb) + 2 * nbi * nbi;
-            } else if ((stag == "HB2ST") || (stag == "SB2ST")) {
+            } else if ((strncmp(stag, "HB2ST", 5) == 0) || (strncmp(stag, "SB2ST", 5) == 0)) {
                 lwork = (2 * nbi + 1) * ni + nbi * nthreads;
             }
-        } else if (algo == "BRD") {
-            if (stag == "2STAG") {
-                lwork = 2 * ni * nbi + ni * max(nbi + 1, factoptnb) + max(2 * nbi * nbi, nbi * nthreads) + (nbi + 1) * ni;
-            } else if (stag == "GE2GB") {
+        } else if (strncmp(algo, "BRD", 3) == 0) {
+            if (strncmp(stag, "2STAG", 5) == 0) {
+                lwork = 2 * ni * nbi + ni * max(nbi + 1, factoptnb) + max((INTEGER)2 * nbi * nbi, nbi * nthreads) + (nbi + 1) * ni;
+            } else if (strncmp(stag, "GE2GB", 5) == 0) {
                 lwork = ni * nbi + ni * max(nbi, factoptnb) + 2 * nbi * nbi;
-            } else if (stag == "GB2BD") {
+            } else if (strncmp(stag, "GB2BD", 5) == 0) {
                 lwork = (3 * nbi + 1) * ni + nbi * nthreads;
             }
         }
