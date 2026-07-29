--- a/mplapack/reference/Cgesvj.cpp
+++ b/mplapack/reference/Cgesvj.cpp
@@ -184,7 +184,11 @@
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
