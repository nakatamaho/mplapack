--- a/mplapack/reference/iMlaenv.cpp
+++ b/mplapack/reference/iMlaenv.cpp
@@ -35,17 +35,20 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
+#include <string.h>
+
+#define subnamlen 17
 
 INTEGER
 iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
     INTEGER return_value = 0;
-    fem::str<16> subnam;
+    char subnam[subnamlen];
+    memset(subnam, '\0', sizeof(subnam));
     INTEGER ic = 0;
-    INTEGER iz = 0;
     INTEGER i = 0;
-    char c1;
     bool sname = false;
     bool cname = false;
+    char c1[1];
     char c2[2];
     char c3[3];
     char c4[2];
@@ -103,63 +106,24 @@
     // Convert NAME to upper case if the first character is lower case.
     //
     return_value = 1;
-    subnam = *name;
-    ic = fem::ichar(subnam(1, 1));
-    iz = fem::ichar("Z");
-    if (iz == 90 || iz == 122) {
-        //
-        // ASCII character set
-        //
-        if (ic >= 97 && ic <= 122) {
-            subnam(1, 1) = fem::fchar(ic - 32);
-            for (i = 2; i <= 6; i = i + 1) {
-                ic = fem::ichar(subnam(i, i));
-                if (ic >= 97 && ic <= 122) {
-                    subnam(i, i) = fem::fchar(ic - 32);
-                }
-            }
-        }
-        //
-    } else if (iz == 233 || iz == 169) {
-        //
-        // EBCDIC character set
-        //
-        if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
-            subnam(1, 1) = fem::fchar(ic + 64);
-            for (i = 2; i <= 6; i = i + 1) {
-                ic = fem::ichar(subnam(i, i));
-                if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
-                    subnam(i, i) = fem::fchar(ic + 64);
-                }
-            }
-        }
-        //
-    } else if (iz == 218 || iz == 250) {
-        //
-        // Prime machines:  ASCII+128
-        //
-        if (ic >= 225 && ic <= 250) {
-            subnam(1, 1) = fem::fchar(ic - 32);
-            for (i = 2; i <= 6; i = i + 1) {
-                ic = fem::ichar(subnam(i, i));
-                if (ic >= 225 && ic <= 250) {
-                    subnam(i, i) = fem::fchar(ic - 32);
-                }
-            }
-        }
-    }
     //
-    c1 = subnam(1, 1);
-    sname = c1 == 'S' || c1 == 'D';
-    cname = c1 == 'C' || c1 == 'Z';
+    // name_len = min((INTEGER)strlen(name), (INTEGER)subnamlen);
+    strcpy(subnam, name);
+    ic = *subnam;
+
+    for (int i = 0; i < strlen(subnam); i++) {
+        subnam[i] = toupper(subnam[i]);
+    }
+    *c1 = *subnam;
+    sname = *c1 == 'R';
+    cname = *c1 == 'C';
     if (!(cname || sname)) {
         return return_value;
     }
-    c2 = subnam(2, 3);
-    c3 = subnam(4, 6);
-    c4[0] = c3[(2 - 1)];
-    c4[1] = c3[(3 - 1)];
-    twostage = fem::len(subnam) >= 11 && strncmp(subnam + 10, "2", 1) == 0;
+    strncpy(c2, subnam + 1, 2);
+    strncpy(c3, subnam + 3, 3);
+    strncpy(c4, c3 + 1, 2);
+    twostage = strlen(subnam) >= 11 && subnam[10] == '2';
     //
     switch (ispec) {
     case 1:
@@ -268,7 +232,7 @@
             } else {
                 nb = 64;
             }
-        } else if (subnam(4, 7) == "QP3RK") {
+        } else if (strncmp(subnam + 3, "QP3RK", 5) == 0) {
             if (sname) {
                 nb = 32;
             } else {
@@ -383,9 +347,9 @@
         } else if (strncmp(c3, "SYL", 3) == 0) {
             // The upper bound is to prevent overly aggressive scaling.
             if (sname) {
-                nb = min(max((INTEGER)48, castINTEGER((min(n1, n2) * 16) / 100)), 240);
+                nb = min(max((INTEGER)48, INTEGER((min(n1, n2) * 16) / 100)), 240);
             } else {
-                nb = min(max((INTEGER)24, castINTEGER((min(n1, n2) * 8) / 100)), 80);
+                nb = min(max((INTEGER)24, INTEGER((min(n1, n2) * 8) / 100)), 80);
             }
         }
     } else if (strncmp(c2, "LA", 2) == 0) {
@@ -449,7 +413,7 @@
             } else {
                 nbmin = 2;
             }
-        } else if (subnam(4, 7) == "QP3RK") {
+        } else if (strncmp(subnam + 3, "QP3RK", 5) == 0) {
             if (sname) {
                 nbmin = 2;
             } else {
@@ -524,7 +488,7 @@
             } else {
                 nx = 128;
             }
-        } else if (subnam(4, 7) == "QP3RK") {
+        } else if (strncmp(subnam + 3, "QP3RK", 5) == 0) {
             if (sname) {
                 nx = 128;
             } else {
