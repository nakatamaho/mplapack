--- a/mplapack/reference/iMparmq.cpp
+++ b/mplapack/reference/iMparmq.cpp
@@ -35,6 +35,9 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <string.h>
+
+#define subnamlen 32
 
 INTEGER
 iMparmq(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER const /* n */, INTEGER const ilo, INTEGER const ihi, INTEGER const /* lwork */) {
@@ -44,6 +47,7 @@
     const INTEGER iacc22 = 16;
     INTEGER nh = 0;
     INTEGER ns = 0;
+    INTEGER name_len;
     const REAL two = 2.0;
     if ((ispec == ishfts) || (ispec == inwin) || (ispec == iacc22)) {
         //
@@ -77,7 +81,7 @@
     const INTEGER inibl = 14;
     const INTEGER nibble = 14;
     const INTEGER knwswp = 500;
-    char subnam[6];
+    char subnam[subnamlen];
     INTEGER ic = 0;
     INTEGER iz = 0;
     INTEGER i = 0;
@@ -129,47 +133,19 @@
         // Convert NAME to upper case if the first character is lower case.
         //
         return_value = 0;
-        subnam = *name;
-        ic = fem::ichar(subnam[0]);
-        iz = fem::ichar("Z");
+        strncpy(subnam, name, subnamlen - 1);
+        ic = *subnam;
+        iz = 'Z';
         if (iz == 90 || iz == 122) {
             //
             // ASCII character set
             //
             if (ic >= 97 && ic <= 122) {
-                subnam[0] = fem::fchar(ic - 32);
-                for (i = 2; i <= 6; i = i + 1) {
-                    ic = fem::ichar(subnam(i, i));
+                *subnam = (char)(ic - 32);
+                for (i = 2; i <= 6; i++) {
+                    ic = subnam[i - 1];
                     if (ic >= 97 && ic <= 122) {
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
-                for (i = 2; i <= 6; i = i + 1) {
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
-                for (i = 2; i <= 6; i = i + 1) {
-                    ic = fem::ichar(subnam(i, i));
-                    if (ic >= 225 && ic <= 250) {
-                        subnam(i, i) = fem::fchar(ic - 32);
+                        subnam[i - 1] = (char)(ic - 32);
                     }
                 }
             }
