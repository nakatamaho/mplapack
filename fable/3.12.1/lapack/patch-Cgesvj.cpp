--- a/mplapack/reference/Cgesvj.cpp
+++ b/mplapack/reference/Cgesvj.cpp
@@ -184,7 +184,11 @@
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
diff --git a/mplapack/reference/Cgesvj.cpp b/mplapack/reference/Cgesvj.cpp
index a701db35..59f4b67a 100644
--- a/mplapack/reference/Cgesvj.cpp
+++ b/mplapack/reference/Cgesvj.cpp
@@ -249,6 +249,12 @@ void Cgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
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
@@ -274,6 +280,12 @@ void Cgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
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
@@ -299,6 +311,12 @@ void Cgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
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
