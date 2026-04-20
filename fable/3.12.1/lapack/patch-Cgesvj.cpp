--- a/mplapack/reference/Cgesvj.cpp
+++ b/mplapack/reference/Cgesvj.cpp
@@ -196,7 +196,11 @@
     sfmin = Rlamch("SafeMinimum");
     rootsfmin = sqrt(sfmin);
     small = sfmin / epsln;
+#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___ ||  defined ___MPLAPACK_BUILD_WITH_MPFR___ ||  defined ___MPLAPACK_BUILD_WITH_GMP___
+    big = one / sfmin;
+#else
     big = Rlamch("Overflow");
+#endif
     // BIG         = ONE    / SFMIN
     rootbig = one / rootsfmin;
     // LARGE = BIG / SQRT( DBLE( M*N ) )
@@ -245,6 +249,12 @@
                 Mxerbla("Cgesvj", -info);
                 return;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            if (aaqq == zero) {
+                sva[p - 1] = zero;
+                continue;
+            }
+#endif
             aaqq = sqrt(aaqq);
             if ((aapp < (big / aaqq)) && noscale) {
                 sva[p - 1] = aapp * aaqq;
@@ -270,6 +280,12 @@
                 Mxerbla("Cgesvj", -info);
                 return;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            if (aaqq == zero) {
+                sva[p - 1] = zero;
+                continue;
+            }
+#endif
             aaqq = sqrt(aaqq);
             if ((aapp < (big / aaqq)) && noscale) {
                 sva[p - 1] = aapp * aaqq;
@@ -295,6 +311,12 @@
                 Mxerbla("Cgesvj", -info);
                 return;
             }
+#if defined ___MPLAPACK_BUILD_WITH_GMP___
+            if (aaqq == zero) {
+                sva[p - 1] = zero;
+                continue;
+            }
+#endif
             aaqq = sqrt(aaqq);
             if ((aapp < (big / aaqq)) && noscale) {
                 sva[p - 1] = aapp * aaqq;
