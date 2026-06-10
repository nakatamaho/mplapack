--- Cpst01.cpp
+++ Cpst01.cpp
@@ -95,7 +95,7 @@
             //
             // Compute the (K,K) element of the result.
             //
-            tr = Cdotc(k, &afac[(k - 1) * ldafac], 1, &afac[(k - 1) * ldafac], 1);
+            tr = Cdotc(k, &afac[(k - 1) * ldafac], 1, &afac[(k - 1) * ldafac], 1).real();
             afac[(k - 1) + (k - 1) * ldafac] = tr;
             //
             // Compute the rest of column K.
