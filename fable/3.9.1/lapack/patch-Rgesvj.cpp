diff --git a/mplapack/reference/Rgesvj.cpp b/mplapack/reference/Rgesvj.cpp
index 493a4877..70758e80 100644
--- a/mplapack/reference/Rgesvj.cpp
+++ b/mplapack/reference/Rgesvj.cpp
@@ -169,7 +169,11 @@ void Rgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
     sfmin = Rlamch("SafeMinimum");
     rootsfmin = sqrt(sfmin);
     small = sfmin / epsln;
+#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
+    big = one / sfmin;
+#else
     big = Rlamch("Overflow");
+#endif
     // BIG         = ONE    / SFMIN
     rootbig = one / rootsfmin;
     large = big / sqrt(castREAL(m * n));
