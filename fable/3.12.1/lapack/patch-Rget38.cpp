--- Rget38.cpp~	2026-04-15 00:00:00.000000000 +0900
+++ Rget38.cpp	2026-04-15 00:00:00.000000000 +0900
@@ -101,6 +101,9 @@ void Rget38(REAL *rmax, INTEGER *lmax, I
     eps = Rlamch("P");
     smlnum = Rlamch("S") / eps;
     bignum = one / smlnum;
+#if defined MPLAPACK_BUILD_WITH_MPFR || defined MPLAPACK_BUILD_WITH_GMP
+    Rlabad(smlnum, bignum);
+#endif
     //
     // EPSIN = 2**(-24) = precision to which input data computed
     //
