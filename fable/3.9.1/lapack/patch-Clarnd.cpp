diff --git a/mplapack/test/matgen/Clarnd.cpp b/mplapack/test/matgen/Clarnd.cpp
index e3d2bebe..b64fd198 100644
--- Clarnd.cpp_	2026-01-18 08:40:08.340948226 +0900
+++ Clarnd.cpp	2026-01-18 08:40:22.211105514 +0900
@@ -51,7 +51,7 @@
     const REAL two = 2.0;
     const REAL one = 1.0;
     const REAL zero = 0.0;
-    const REAL twopi = 6.28318530717958647692528676655900576839;
+    const REAL twopi = two * pi(zero);
     if (idist == 1) {
         //
         // real and imaginary parts each uniform (0,1)
