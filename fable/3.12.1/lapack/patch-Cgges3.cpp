diff --git a/mplapack/reference/Cgges3.cpp b/mplapack/reference/Cgges3.cpp
index 71a7f4f80..e971ad96b 100644
--- a/mplapack/reference/Cgges3.cpp
+++ b/mplapack/reference/Cgges3.cpp
@@ -173,6 +173,9 @@ void Cgges3(const char *jobvsl, const char *jobvsr, const char *sort, bool (*sel
     eps = Rlamch("P");
     smlnum = Rlamch("S");
     bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_MPFR___ ||  defined ___MPLAPACK_BUILD_WITH_GMP___ ||  defined ___MPLAPACK_BUILD_WITH_BINARY80___ ||  defined ___MPLAPACK_BUILD_WITH_BINARY128___
+    Rlabad(smlnum, bignum);
+#endif
     smlnum = sqrt(smlnum) / eps;
     bignum = one / smlnum;
     //
