--- Cchkeq.cpp
+++ Cchkeq.cpp
@@ -60,12 +60,12 @@
     const INTEGER nsz = 5;
     const INTEGER npow = 2 * nsz + 1;
     const REAL ten = 10.0;
-    REAL pow[npow];
+    REAL pow10[npow];
     const REAL one = 1.0;
     REAL rpow[npow];
     for (i = 1; i <= npow; i = i + 1) {
-        pow[i - 1] = pow(ten, (i - 1));
-        rpow[i - 1] = one / pow[i - 1];
+        pow10[i - 1] = pow(ten, (i - 1));
+        rpow[i - 1] = one / pow10[i - 1];
     }
     //
     // Test Cgeequ
@@ -87,7 +87,7 @@
             for (j = 1; j <= nsz; j = j + 1) {
                 for (i = 1; i <= nsz; i = i + 1) {
                     if (i <= m && j <= n) {
-                        a[(i - 1) + (j - 1) * nsz] = pow[(i + j + 1) - 1] * pow((-1), (i + j));
+                        a[(i - 1) + (j - 1) * nsz] = pow10[(i + j + 1) - 1] * ((((i + j) % 2) == 0) ? one : -one);
                     } else {
                         a[(i - 1) + (j - 1) * nsz] = czero;
                     }
@@ -102,12 +102,12 @@
                 if (n != 0 && m != 0) {
                     reslts[1 - 1] = max(reslts[1 - 1], abs((rcond - rpow[m - 1]) / rpow[m - 1]));
                     reslts[1 - 1] = max(reslts[1 - 1], abs((ccond - rpow[n - 1]) / rpow[n - 1]));
-                    reslts[1 - 1] = max(reslts[1 - 1], abs((norm - pow[(n + m + 1) - 1]) / pow[(n + m + 1) - 1]));
+                    reslts[1 - 1] = max(reslts[1 - 1], abs((norm - pow10[(n + m + 1) - 1]) / pow10[(n + m + 1) - 1]));
                     for (i = 1; i <= m; i = i + 1) {
                         reslts[1 - 1] = max(reslts[1 - 1], abs((r[i - 1] - rpow[(i + n + 1) - 1]) / rpow[(i + n + 1) - 1]));
                     }
                     for (j = 1; j <= n; j = j + 1) {
-                        reslts[1 - 1] = max(reslts[1 - 1], abs((c[j - 1] - pow[(n - j + 1) - 1]) / pow[(n - j + 1) - 1]));
+                        reslts[1 - 1] = max(reslts[1 - 1], abs((c[j - 1] - pow10[(n - j + 1) - 1]) / pow10[(n - j + 1) - 1]));
                     }
                 }
             }
@@ -160,7 +160,7 @@
                     for (j = 1; j <= n; j = j + 1) {
                         for (i = 1; i <= m; i = i + 1) {
                             if (i <= min(m, j + kl) && i >= max((INTEGER)1, j - ku) && j <= n) {
-                                ab[((ku + 1 + i - j) - 1) + (j - 1) * nszb] = pow[(i + j + 1) - 1] * pow((-1), (i + j));
+                                ab[((ku + 1 + i - j) - 1) + (j - 1) * nszb] = pow10[(i + j + 1) - 1] * ((((i + j) % 2) == 0) ? one : -one);
                             }
                         }
                     }
@@ -192,12 +192,12 @@
                             ratio = rcmin / rcmax;
                             reslts[2 - 1] = max(reslts[2 - 1], abs((ccond - ratio) / ratio));
                             //
-                            reslts[2 - 1] = max(reslts[2 - 1], abs((norm - pow[(n + m + 1) - 1]) / pow[(n + m + 1) - 1]));
+                            reslts[2 - 1] = max(reslts[2 - 1], abs((norm - pow10[(n + m + 1) - 1]) / pow10[(n + m + 1) - 1]));
                             for (i = 1; i <= m; i = i + 1) {
                                 rcmax = zero;
                                 for (j = 1; j <= n; j = j + 1) {
                                     if (i <= j + kl && i >= j - ku) {
-                                        ratio = abs(r[i - 1] * pow[(i + j + 1) - 1] * c[j - 1]);
+                                        ratio = abs(r[i - 1] * pow10[(i + j + 1) - 1] * c[j - 1]);
                                         rcmax = max(rcmax, ratio);
                                     }
                                 }
@@ -208,7 +208,7 @@
                                 rcmax = zero;
                                 for (i = 1; i <= m; i = i + 1) {
                                     if (i <= j + kl && i >= j - ku) {
-                                        ratio = abs(r[i - 1] * pow[(i + j + 1) - 1] * c[j - 1]);
+                                        ratio = abs(r[i - 1] * pow10[(i + j + 1) - 1] * c[j - 1]);
                                         rcmax = max(rcmax, ratio);
                                     }
                                 }
@@ -230,7 +230,7 @@
         for (i = 1; i <= nsz; i = i + 1) {
             for (j = 1; j <= nsz; j = j + 1) {
                 if (i <= n && j == i) {
-                    a[(i - 1) + (j - 1) * nsz] = pow[(i + j + 1) - 1] * pow((-1), (i + j));
+                    a[(i - 1) + (j - 1) * nsz] = pow10[(i + j + 1) - 1] * ((((i + j) % 2) == 0) ? one : -one);
                 } else {
                     a[(i - 1) + (j - 1) * nsz] = czero;
                 }
@@ -244,7 +244,7 @@
         } else {
             if (n != 0) {
                 reslts[3 - 1] = max(reslts[3 - 1], abs((rcond - rpow[n - 1]) / rpow[n - 1]));
-                reslts[3 - 1] = max(reslts[3 - 1], abs((norm - pow[(2 * n + 1) - 1]) / pow[(2 * n + 1) - 1]));
+                reslts[3 - 1] = max(reslts[3 - 1], abs((norm - pow10[(2 * n + 1) - 1]) / pow10[(2 * n + 1) - 1]));
                 for (i = 1; i <= n; i = i + 1) {
                     reslts[3 - 1] = max(reslts[3 - 1], abs((r[i - 1] - rpow[(i + 1) - 1]) / rpow[(i + 1) - 1]));
                 }
@@ -270,7 +270,7 @@
             ap[i - 1] = czero;
         }
         for (i = 1; i <= n; i = i + 1) {
-            ap[((i * (i + 1)) / 2) - 1] = pow[(2 * i + 1) - 1];
+            ap[((i * (i + 1)) / 2) - 1] = pow10[(2 * i + 1) - 1];
         }
         //
         Cppequ("U", n, ap, r, rcond, norm, info);
@@ -280,7 +280,7 @@
         } else {
             if (n != 0) {
                 reslts[4 - 1] = max(reslts[4 - 1], abs((rcond - rpow[n - 1]) / rpow[n - 1]));
-                reslts[4 - 1] = max(reslts[4 - 1], abs((norm - pow[(2 * n + 1) - 1]) / pow[(2 * n + 1) - 1]));
+                reslts[4 - 1] = max(reslts[4 - 1], abs((norm - pow10[(2 * n + 1) - 1]) / pow10[(2 * n + 1) - 1]));
                 for (i = 1; i <= n; i = i + 1) {
                     reslts[4 - 1] = max(reslts[4 - 1], abs((r[i - 1] - rpow[(i + 1) - 1]) / rpow[(i + 1) - 1]));
                 }
@@ -294,7 +294,7 @@
         }
         j = 1;
         for (i = 1; i <= n; i = i + 1) {
-            ap[j - 1] = pow[(2 * i + 1) - 1];
+            ap[j - 1] = pow10[(2 * i + 1) - 1];
             j += (n - i + 1);
         }
         //
@@ -305,7 +305,7 @@
         } else {
             if (n != 0) {
                 reslts[4 - 1] = max(reslts[4 - 1], abs((rcond - rpow[n - 1]) / rpow[n - 1]));
-                reslts[4 - 1] = max(reslts[4 - 1], abs((norm - pow[(2 * n + 1) - 1]) / pow[(2 * n + 1) - 1]));
+                reslts[4 - 1] = max(reslts[4 - 1], abs((norm - pow10[(2 * n + 1) - 1]) / pow10[(2 * n + 1) - 1]));
                 for (i = 1; i <= n; i = i + 1) {
                     reslts[4 - 1] = max(reslts[4 - 1], abs((r[i - 1] - rpow[(i + 1) - 1]) / rpow[(i + 1) - 1]));
                 }
@@ -334,7 +334,7 @@
                 }
             }
             for (j = 1; j <= n; j = j + 1) {
-                ab[((kl + 1) - 1) + (j - 1) * nszb] = pow[(2 * j + 1) - 1];
+                ab[((kl + 1) - 1) + (j - 1) * nszb] = pow10[(2 * j + 1) - 1];
             }
             //
             Cpbequ("U", n, kl, ab, nszb, r, rcond, norm, info);
@@ -344,7 +344,7 @@
             } else {
                 if (n != 0) {
                     reslts[5 - 1] = max(reslts[5 - 1], abs((rcond - rpow[n - 1]) / rpow[n - 1]));
-                    reslts[5 - 1] = max(reslts[5 - 1], abs((norm - pow[(2 * n + 1) - 1]) / pow[(2 * n + 1) - 1]));
+                    reslts[5 - 1] = max(reslts[5 - 1], abs((norm - pow10[(2 * n + 1) - 1]) / pow10[(2 * n + 1) - 1]));
                     for (i = 1; i <= n; i = i + 1) {
                         reslts[5 - 1] = max(reslts[5 - 1], abs((r[i - 1] - rpow[(i + 1) - 1]) / rpow[(i + 1) - 1]));
                     }
@@ -366,7 +366,7 @@
                 }
             }
             for (j = 1; j <= n; j = j + 1) {
-                ab[(j - 1) * nszb] = pow[(2 * j + 1) - 1];
+                ab[(j - 1) * nszb] = pow10[(2 * j + 1) - 1];
             }
             //
             Cpbequ("L", n, kl, ab, nszb, r, rcond, norm, info);
@@ -376,7 +376,7 @@
             } else {
                 if (n != 0) {
                     reslts[5 - 1] = max(reslts[5 - 1], abs((rcond - rpow[n - 1]) / rpow[n - 1]));
-                    reslts[5 - 1] = max(reslts[5 - 1], abs((norm - pow[(2 * n + 1) - 1]) / pow[(2 * n + 1) - 1]));
+                    reslts[5 - 1] = max(reslts[5 - 1], abs((norm - pow10[(2 * n + 1) - 1]) / pow10[(2 * n + 1) - 1]));
                     for (i = 1; i <= n; i = i + 1) {
                         reslts[5 - 1] = max(reslts[5 - 1], abs((r[i - 1] - rpow[(i + 1) - 1]) / rpow[(i + 1) - 1]));
                     }
