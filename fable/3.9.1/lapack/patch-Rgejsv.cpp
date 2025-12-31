diff --git a/mplapack/reference/Rgejsv.cpp b/mplapack/reference/Rgejsv.cpp
index 89dca753..1e95a19d 100644
--- a/mplapack/reference/Rgejsv.cpp
+++ b/mplapack/reference/Rgejsv.cpp
@@ -154,7 +154,11 @@ void Rgejsv(const char *joba, const char *jobu, const char *jobv, const char *jo
     epsln = Rlamch("Epsilon");
     sfmin = Rlamch("SafeMinimum");
     small = sfmin / epsln;
+#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
+    big = one / sfmin;
+#else
     big = Rlamch("O");
+#endif
     // BIG   = ONE / SFMIN
     //
     // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
