diff --git a/mplapack/reference/Cgges.cpp b/mplapack/reference/Cgges.cpp
index 000000000..000000000 100644
--- a/mplapack/reference/Cgges.cpp
+++ b/mplapack/reference/Cgges.cpp
@@ -164,6 +164,9 @@ void Cgges(const char *jobvsl, const char *jobvsr, const char *sort, bool (*selc
     eps = Rlamch("P");
     smlnum = Rlamch("S");
     bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
+    Rlabad(smlnum, bignum);
+#endif
     smlnum = sqrt(smlnum) / eps;
     bignum = one / smlnum;
     //
