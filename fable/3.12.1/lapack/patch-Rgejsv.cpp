--- a/mplapack/reference/Rgejsv.cpp
+++ b/mplapack/reference/Rgejsv.cpp
@@ -161,7 +161,11 @@
     epsln = Rlamch("Epsilon");
     sfmin = Rlamch("SafeMinimum");
     small = sfmin / epsln;
+#if defined MPLAPACK_BUILD_WITH_DD || defined MPLAPACK_BUILD_WITH_QD ||  defined MPLAPACK_BUILD_WITH_MPFR ||  defined MPLAPACK_BUILD_WITH_GMP
+    big = one / sfmin;
+#else
     big = Rlamch("O");
+#endif
     // BIG   = ONE / SFMIN
     //
     // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
@@ -182,6 +186,12 @@
             Mxerbla("Rgejsv", -info);
             return;
         }
+#if defined MPLAPACK_BUILD_WITH_GMP
+        if (aaqq == zero) {
+            sva[p - 1] = zero;
+            continue;
+        }
+#endif
         aaqq = sqrt(aaqq);
         if ((aapp < (big / aaqq)) && noscal) {
             sva[p - 1] = aapp * aaqq;
