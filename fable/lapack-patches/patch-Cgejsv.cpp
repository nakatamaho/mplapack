diff --git a/mplapack/reference/Cgejsv.cpp b/mplapack/reference/Cgejsv.cpp
index e1dc6043..74575006 100644
--- a/mplapack/reference/Cgejsv.cpp
+++ b/mplapack/reference/Cgejsv.cpp
@@ -409,11 +409,7 @@ void Cgejsv(const char *joba, const char *jobu, const char *jobv, const char *jo
     epsln = Rlamch("Epsilon");
     sfmin = Rlamch("SafeMinimum");
     small = sfmin / epsln;
-#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
-    big = one / sfmin;
-#else
-    big = Rlamch("Overflow");
-#endif
+    big = Rlamch("O");
     // BIG   = ONE / SFMIN
     //
     // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
