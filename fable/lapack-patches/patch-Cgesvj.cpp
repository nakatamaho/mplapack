diff --git a/mplapack/reference/Cgesvj.cpp b/mplapack/reference/Cgesvj.cpp
index bf60a606..78bf270e 100644
--- a/mplapack/reference/Cgesvj.cpp
+++ b/mplapack/reference/Cgesvj.cpp
@@ -177,11 +179,7 @@ void Cgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const
     sfmin = Rlamch("SafeMinimum");
     rootsfmin = sqrt(sfmin);
     small = sfmin / epsln;
-#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
-    big = one / sfmin;
-#else
     big = Rlamch("Overflow");
-#endif
     // BIG         = ONE    / SFMIN
     rootbig = one / rootsfmin;
     // LARGE = BIG / SQRT( DBLE( M*N ) )
