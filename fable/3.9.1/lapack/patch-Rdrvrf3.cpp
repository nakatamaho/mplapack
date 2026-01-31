--- Rdrvrf3.cpp
+++ Rdrvrf3.cpp
@@ -80,9 +80,9 @@
     INTEGER idiag = 0;
     fem::str<1> diag;
     INTEGER ialpha = 0;
-    const REAL zero = COMPLEX(0.0, 0.0);
+    const REAL zero = 0.0;
     REAL alpha = 0.0;
-    const REAL one = COMPLEX(1.0, 0.0);
+    const REAL one = 1.0;
     INTEGER na = 0;
     INTEGER j = 0;
     const INTEGER ntests = 1;
