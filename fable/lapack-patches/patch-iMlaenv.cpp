diff --git a/mplapack/reference/iMlaenv.cpp b/mplapack/reference/iMlaenv.cpp
index f363fa01..5ff67569 100644
--- a/mplapack/reference/iMlaenv.cpp
+++ b/mplapack/reference/iMlaenv.cpp
@@ -28,20 +28,17 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
-#include <string.h>
-
-#define subnamlen 17
 
 INTEGER
 iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
     INTEGER return_value = 0;
-    char subnam[subnamlen];
-    memset(subnam, '\0', sizeof(subnam));
+    fem::str<16> subnam = 0;
     INTEGER ic = 0;
+    INTEGER iz = 0;
     INTEGER i = 0;
+    char c1;
     bool sname = false;
     bool cname = false;
-    char c1[1];
     char c2[2];
     char c3[3];
     char c4[2];
@@ -49,7 +46,6 @@ iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n
     INTEGER nb = 0;
     INTEGER nbmin = 0;
     INTEGER nx = 0;
-    INTEGER name_len;
     //
     switch (ispec) {
     case 1:
@@ -98,24 +94,62 @@ statement_10:
     // Convert NAME to upper case if the first character is lower case.
     //
     return_value = 1;
-    //
-    // name_len = min((INTEGER)strlen(name), (INTEGER)subnamlen);
-    strcpy(subnam, name);
-    ic = *subnam;
-
-    for (int i = 0; i < strlen(subnam); i++) {
-        subnam[i] = toupper(subnam[i]);
+    subnam = name;
+    ic = fem::ichar(subnam(1, 1));
+    iz = fem::ichar("Z");
+    if (iz == 90 || iz == 122) {
+        //
+        // ASCII character set
+        //
+        if (ic >= 97 && ic <= 122) {
+            subnam(1, 1) = fem::fchar(ic - 32);
+            for (i = 2; i <= 6; i = i + 1) {
+                ic = fem::ichar(subnam(i, i));
+                if (ic >= 97 && ic <= 122) {
+                    subnam(i, i) = fem::fchar(ic - 32);
+                }
+            }
+        }
+        //
+    } else if (iz == 233 || iz == 169) {
+        //
+        // EBCDIC character set
+        //
+        if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
+            subnam(1, 1) = fem::fchar(ic + 64);
+            for (i = 2; i <= 6; i = i + 1) {
+                ic = fem::ichar(subnam(i, i));
+                if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
+                    subnam(i, i) = fem::fchar(ic + 64);
+                }
+            }
+        }
+        //
+    } else if (iz == 218 || iz == 250) {
+        //
+        // Prime machines:  ASCII+128
+        //
+        if (ic >= 225 && ic <= 250) {
+            subnam(1, 1) = fem::fchar(ic - 32);
+            for (i = 2; i <= 6; i = i + 1) {
+                ic = fem::ichar(subnam(i, i));
+                if (ic >= 225 && ic <= 250) {
+                    subnam(i, i) = fem::fchar(ic - 32);
+                }
+            }
+        }
     }
-    *c1 = *subnam;
-    sname = *c1 == 'R';
-    cname = *c1 == 'C';
+    //
+    c1 = subnam(1, 1);
+    sname = c1 == 'S' || c1 == 'D';
+    cname = c1 == 'C' || c1 == 'Z';
     if (!(cname || sname)) {
         return return_value;
     }
-    strncpy(c2, subnam + 1, 2);
-    strncpy(c3, subnam + 3, 3);
-    strncpy(c4, c3 + 1, 2);
-    twostage = strlen(subnam) >= 11 && subnam[10] == '2';
+    c2 = subnam(2, 3);
+    c3 = subnam(4, 6);
+    c4 = c3(2, 3);
+    twostage = fem::len(subnam) >= 11 && strncmp(subnam + 10, "2", 1) == 0;
     //
     switch (ispec) {
     case 1:
