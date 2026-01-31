--- Cpbt01.cpp
+++ Cpbt01.cpp
@@ -96,7 +96,7 @@
             //
             // Compute the (K,K) element of the result.
             //
-            akk = Cdotc(klen + 1, &afac[(kc - 1) + (k - 1) * ldafac], 1, &afac[(kc - 1) + (k - 1) * ldafac], 1);
+            akk = Cdotc(klen + 1, &afac[(kc - 1) + (k - 1) * ldafac], 1, &afac[(kc - 1) + (k - 1) * ldafac], 1).real();
             afac[((kd + 1) - 1) + (k - 1) * ldafac] = akk;
             //
             // Compute the rest of column K.
