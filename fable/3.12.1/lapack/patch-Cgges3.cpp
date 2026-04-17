--- a/mplapack/reference/Cgges3.cpp
+++ b/mplapack/reference/Cgges3.cpp
@@ -173,6 +173,9 @@
     eps = Rlamch("P");
     smlnum = Rlamch("S");
     bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___ ||  defined ___MPLAPACK_BUILD_WITH_GMP___ ||  defined ___MPLAPACK_BUILD_WITH_BINARY80___ ||  defined ___MPLAPACK_BUILD_WITH_BINARY128___
+    Rlabad(smlnum, bignum);
+#endif
     smlnum = sqrt(smlnum) / eps;
     bignum = one / smlnum;
     //
