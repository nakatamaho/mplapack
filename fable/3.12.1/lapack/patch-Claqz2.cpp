diff --git a/mplapack/reference/Claqz2.cpp b/mplapack/reference/Claqz2.cpp
--- a/mplapack/reference/Claqz2.cpp
+++ b/mplapack/reference/Claqz2.cpp
@@ -168,8 +168,10 @@ void Claqz2(bool const ilschur, bool const ilq, bool const ilz, INTEGER const n,
     INTEGER istopm = 0;
     if (kwtop != ilo && s != czero) {
         // Reflect spike back, this will create optimally packed bulges
-        // FABLE_UNHANDLED_SLICE_ASSIGN: a[__SLICE2D__(kwtop, kwbot, kwtop - 1, lda)] = a[(kwtop - 1) + ((kwtop - 1) - 1)*lda] * conj(qc[__FABLE_UNHANDLED_SLICE2D_COL__(1, 1, jw - nd, ldqc)]);
-        FABLE_UNHANDLED_SLICE_ASSIGN();
+        const COMPLEX spike = a[(kwtop - 1) + ((kwtop - 1) - 1) * lda];
+        for (INTEGER i_ = kwtop; i_ <= kwbot; i_++) {
+            a[(i_ - 1) + ((kwtop - 1) - 1) * lda] = spike * conj(qc[(i_ - kwtop) * ldqc]);
+        }
         for (k = kwbot - 1; k >= kwtop; k = k - 1) {
             Clartg(a[(k - 1) + ((kwtop - 1) - 1) * lda], a[((k + 1) - 1) + ((kwtop - 1) - 1) * lda], c1, s1, temp);
             a[(k - 1) + ((kwtop - 1) - 1) * lda] = temp;
