--- a/mplapack/reference/Cgeev.cpp
+++ b/mplapack/reference/Cgeev.cpp
@@ -149,6 +149,9 @@
     eps = Rlamch("P");
     smlnum = Rlamch("S");
     bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
+    Rlabad(smlnum, bignum);
+#endif
     smlnum = sqrt(smlnum) / eps;
     bignum = one / smlnum;
     //
