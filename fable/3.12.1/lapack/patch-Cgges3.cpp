--- a/mplapack/reference/Cgges3.cpp
+++ b/mplapack/reference/Cgges3.cpp
@@ -173,6 +173,9 @@
     eps = Rlamch("P");
     smlnum = Rlamch("S");
     bignum = one / smlnum;
+#if defined MPLAPACK_BUILD_WITH_MPFR ||  defined MPLAPACK_BUILD_WITH_GMP ||  defined MPLAPACK_BUILD_WITH_BINARY80 ||  defined MPLAPACK_BUILD_WITH_BINARY128
+    Rlabad(smlnum, bignum);
+#endif
     smlnum = sqrt(smlnum) / eps;
     bignum = one / smlnum;
     //
