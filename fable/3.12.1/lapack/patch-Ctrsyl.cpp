--- a/mplapack/reference/Ctrsyl.cpp
+++ b/mplapack/reference/Ctrsyl.cpp
@@ -84,7 +84,7 @@
     bignum = one / smlnum;
     REAL dum[1];
     REAL smin = max(smlnum, eps * Clange("M", m, m, a, lda, dum), eps * Clange("M", n, n, b, ldb, dum));
-    REAL sgn = isgn;
+    REAL sgn = castREAL(isgn);
     //
     INTEGER l = 0;
     INTEGER k = 0;
