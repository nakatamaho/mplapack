--- Cerrlqt.cpp
+++ Cerrlqt.cpp
@@ -55,15 +55,16 @@
     INTEGER j = 0;
     const INTEGER nmax = 2;
     INTEGER i = 0;
+    const REAL one = 1.0;
     COMPLEX a[nmax * nmax];
     COMPLEX c[nmax * nmax];
     COMPLEX t[nmax * nmax];
     COMPLEX w[nmax];
     for (j = 1; j <= nmax; j = j + 1) {
         for (i = 1; i <= nmax; i = i + 1) {
-            a[(i - 1) + (j - 1) * nmax] = 1.0 / COMPLEX(castREAL(i + j), 0.0);
-            c[(i - 1) + (j - 1) * nmax] = 1.0 / COMPLEX(castREAL(i + j), 0.0);
-            t[(i - 1) + (j - 1) * nmax] = 1.0 / COMPLEX(castREAL(i + j), 0.0);
+            a[(i - 1) + (j - 1) * nmax] = one / COMPLEX(castREAL(i + j), 0.0);
+            c[(i - 1) + (j - 1) * nmax] = one / COMPLEX(castREAL(i + j), 0.0);
+            t[(i - 1) + (j - 1) * nmax] = one / COMPLEX(castREAL(i + j), 0.0);
         }
         w[j - 1] = 0.0;
     }
