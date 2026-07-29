--- a/mplapack/reference/Rgesvj.cpp
+++ b/mplapack/reference/Rgesvj.cpp
@@ -188,7 +188,11 @@
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
     large = big / sqrt(castREAL(m * n));
@@ -237,6 +241,12 @@
                 Mxerbla("Rgesvj", -info);
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
@@ -262,6 +272,12 @@
                 Mxerbla("Rgesvj", -info);
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
@@ -287,6 +303,12 @@
                 Mxerbla("Rgesvj", -info);
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
