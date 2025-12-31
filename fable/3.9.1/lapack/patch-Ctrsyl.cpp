diff --git a/mplapack/reference/Ctrsyl.cpp b/mplapack/reference/Ctrsyl.cpp
index c5e71eeb..d9e982eb 100644
--- a/mplapack/reference/Ctrsyl.cpp
+++ b/mplapack/reference/Ctrsyl.cpp
@@ -84,7 +84,7 @@ void Ctrsyl(const char *trana, const char *tranb, INTEGER const isgn, INTEGER co
     bignum = one / smlnum;
     REAL dum[1];
     REAL smin = max(smlnum, eps * Clange("M", m, m, a, lda, dum), eps * Clange("M", n, n, b, ldb, dum));
-    REAL sgn = isgn;
+    REAL sgn = castREAL(isgn);
     //
     INTEGER l = 0;
     INTEGER k = 0;
