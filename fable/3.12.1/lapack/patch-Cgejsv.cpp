--- a/mplapack/reference/Cgejsv.cpp
+++ b/mplapack/reference/Cgejsv.cpp
@@ -415,7 +415,11 @@
     epsln = Rlamch("Epsilon");
     sfmin = Rlamch("SafeMinimum");
     small = sfmin / epsln;
+#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___ ||  defined ___MPLAPACK_BUILD_WITH_MPFR___ ||  defined ___MPLAPACK_BUILD_WITH_GMP___
+    big = one / sfmin;
+#else
     big = Rlamch("O");
+#endif
     // BIG   = ONE / SFMIN
     //
     // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
@@ -436,6 +440,12 @@
             Mxerbla("Cgejsv", -info);
             return;
         }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+        if (aaqq == zero) {
+            sva[p - 1] = zero;
+            continue;
+        }
+#endif
         aaqq = sqrt(aaqq);
         if ((aapp < (big / aaqq)) && noscal) {
             sva[p - 1] = aapp * aaqq;
