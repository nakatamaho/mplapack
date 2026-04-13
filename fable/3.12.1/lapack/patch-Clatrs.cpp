--- a/mplapack/reference/Clatrs.cpp
+++ b/mplapack/reference/Clatrs.cpp
@@ -36,6 +36,8 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+inline REAL cabs2(COMPLEX zdum) { return abs(zdum.real() / 2.0) + abs(zdum.imag() / 2.0); }
+
 void Clatrs(const char *uplo, const char *trans, const char *diag, const char *normin, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, REAL &scale, REAL *cnorm, INTEGER &info) {
     COMPLEX zdum = 0.0;
     bool upper = false;

diff --git a/mplapack/reference/Clatrs.cpp b/mplapack/reference/Clatrs.cpp
index f7608d507..0d366183b 100644
--- a/mplapack/reference/Clatrs.cpp
+++ b/mplapack/reference/Clatrs.cpp
@@ -94,15 +94,20 @@ void Clatrs(const char *uplo, const char *trans, const char *diag, const char *n
     //
     // Quick return if possible
     //
-    scale = one;
     if (n == 0) {
         return;
     }
     //
     // Determine machine dependent parameters to control overflow.
     //
-    smlnum = Rlamch("Safe minimum") / Rlamch("Precision");
+    smlnum = Rlamch("Safe minimum");
+    bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_BINARY128___ ||  defined ___MPLAPACK_BUILD_WITH_BINARY80___
+    Rlabad(smlnum, bignum);
+#endif
+    smlnum = smlnum / Rlamch("Precision");
     bignum = one / smlnum;
+    scale = one;
     //
     if (Mlsame(normin, "N")) {
         //
