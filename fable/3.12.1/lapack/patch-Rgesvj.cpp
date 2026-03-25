--- a/mplapack/reference/Rgesvj.cpp
+++ b/mplapack/reference/Rgesvj.cpp
@@ -174,7 +174,11 @@
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
     large = big / sqrt(castREAL(m * n));
diff --git a/mplapack/reference/Rgesvj.cpp b/mplapack/reference/Rgesvj.cpp
index 7415c975..967768f2 100644
--- a/mplapack/reference/Rgesvj.cpp
+++ b/mplapack/reference/Rgesvj.cpp
@@ -241,6 +241,12 @@ void Rgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
                 Mxerbla("Rgesvj", -info);
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
@@ -266,6 +272,12 @@ void Rgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
                 Mxerbla("Rgesvj", -info);
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
@@ -291,6 +303,12 @@ void Rgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
                 Mxerbla("Rgesvj", -info);
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
