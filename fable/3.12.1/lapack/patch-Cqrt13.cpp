--- Cqrt13.cpp	2026-04-12 18:51:10.001470517 +0900
+++ Cqrt13.cpp	2026-04-12 18:54:48.950639775 +0900
@@ -70,6 +70,9 @@
         norma = Clange("Max", m, n, a, lda, dummy);
         smlnum = Rlamch("Safe minimum");
         bignum = one / smlnum;
+#if defined MPLAPACK_BUILD_WITH_BINARY80 || defined MPLAPACK_BUILD_WITH_BINARY128
+        Rlabad(smlnum, bignum);
+#endif
         smlnum = smlnum / Rlamch("Epsilon");
         bignum = one / smlnum;
         //
