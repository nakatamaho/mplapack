diff --git a/mplapack/test/lin/common/Cqrt13.cpp b/mplapack/test/lin/common/Cqrt13.cpp
index 658d39e3..acc879af 100644
--- Cqrt13.cpp
+++ Cqrt13.cpp
@@ -72,6 +72,20 @@ void Cqrt13(INTEGER const scale, INTEGER const m, INTEGER const n, COMPLEX *a, I
         bignum = one / smlnum;
         smlnum = smlnum / Rlamch("Epsilon");
         bignum = one / smlnum;
+#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
+        // Historical DLABAD-style range reduction for very wide exponent ranges.
+        //
+        // Current LAPACK DLABAD is a no-op, but its former logic reduced SMALL and
+        // LARGE by square roots when LOG10(LARGE) > 2000. That condition is false
+        // for ordinary IEEE binary32/binary64 and true for wide-range formats such
+        // as binary80/binary128. Apply the same policy here so that standard
+        // backends keep LAPACK's original test scaling while wide-range backends
+        // avoid excessively extreme near-underflow/near-overflow matrices.
+        if (log10(bignum) > 2000.0) {
+            smlnum = sqrt(smlnum);
+            bignum = sqrt(bignum);
+        }
+#endif
         //
         if (scale == 2) {
             //
