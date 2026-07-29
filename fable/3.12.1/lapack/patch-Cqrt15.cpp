--- Cqrt15.cpp	2026-04-12 18:40:34.639159027 +0900
+++ Cqrt15.cpp	2026-04-12 18:41:33.512871309 +0900
@@ -67,6 +67,9 @@
     //
     smlnum = Rlamch("Safe minimum");
     bignum = one / smlnum;
+#if defined MPLAPACK_BUILD_WITH_BINARY80 || defined MPLAPACK_BUILD_WITH_BINARY128
+    Rlabad(smlnum, bignum);
+#endif
     eps = Rlamch("Epsilon");
     smlnum = (smlnum / eps) / eps;
     bignum = one / smlnum;
