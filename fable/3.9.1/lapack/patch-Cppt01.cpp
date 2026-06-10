--- Cppt01.cpp
+++ Cppt01.cpp
@@ -97,7 +97,7 @@
             //
             // Compute the (K,K) element of the result.
             //
-            tr = Cdotc(k, &afac[kc - 1], 1, &afac[kc - 1], 1);
+	    tr = Cdotc(k, &afac[kc - 1], 1, &afac[kc - 1], 1).real();
             afac[(kc + k - 1) - 1] = tr;
             //
             // Compute the rest of column K.
