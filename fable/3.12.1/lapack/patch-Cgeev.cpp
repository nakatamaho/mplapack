diff --git a/mplapack/reference/Cgeev.cpp b/mplapack/reference/Cgeev.cpp
index 000000000..000000000 100644
--- a/mplapack/reference/Cgeev.cpp
+++ b/mplapack/reference/Cgeev.cpp
@@ -149,6 +149,9 @@ void Cgeev(const char *jobvl, const char *jobvr, INTEGER const n, COMPLEX *a, IN
     eps = Rlamch("P");
     smlnum = Rlamch("S");
     bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
+    Rlabad(smlnum, bignum);
+#endif
     smlnum = sqrt(smlnum) / eps;
     bignum = one / smlnum;
     //
