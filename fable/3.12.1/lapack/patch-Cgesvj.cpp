--- a/mplapack/reference/Cgesvj.cpp
+++ b/mplapack/reference/Cgesvj.cpp
@@ -196,7 +196,11 @@
     sfmin = Rlamch("SafeMinimum");
     rootsfmin = sqrt(sfmin);
     small = sfmin / epsln;
+#if defined MPLAPACK_BUILD_WITH_DD || defined MPLAPACK_BUILD_WITH_QD ||  defined MPLAPACK_BUILD_WITH_MPFR ||  defined MPLAPACK_BUILD_WITH_GMP
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
+#if defined MPLAPACK_BUILD_WITH_GMP
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
+#if defined MPLAPACK_BUILD_WITH_GMP
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
+#if defined MPLAPACK_BUILD_WITH_GMP
+            if (aaqq == zero) {
+                sva[p - 1] = zero;
+                continue;
+            }
+#endif
             aaqq = sqrt(aaqq);
             if ((aapp < (big / aaqq)) && noscale) {
                 sva[p - 1] = aapp * aaqq;
