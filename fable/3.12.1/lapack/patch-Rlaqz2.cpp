--- a/mplapack/reference/Rlaqz2.cpp
+++ b/mplapack/reference/Rlaqz2.cpp
@@ -35,15 +35,6 @@
 
 void Rlaqz2(bool const ilq, bool const ilz, INTEGER const k, INTEGER const istartm, INTEGER const istopm, INTEGER const ihi, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER const nq, INTEGER const qstart, REAL *q, INTEGER const ldq, INTEGER const nz, INTEGER const zstart, REAL *z, INTEGER const ldz) {
     INTEGER ldh = 2;
-    //
-    // Arguments
-    //
-    // Parameters
-    //
-    // Local variables
-    //
-    // External functions
-    //
     REAL h[2 * 3];
     REAL c1 = 0.0;
     REAL s1 = 0.0;
@@ -53,7 +44,13 @@
     REAL s2 = 0.0;
     if (k + 2 == ihi) {
         // Shift is located on the edge of the matrix, remove it
-        h = b(ihi - 1, ihi, ihi - 2, ihi);
+        // H = B( IHI-1:IHI, IHI-2:IHI )  -- 2x3 submatrix
+        h[(1 - 1) + (1 - 1) * 2] = b[((ihi - 1) - 1) + ((ihi - 2) - 1) * ldb];
+        h[(2 - 1) + (1 - 1) * 2] = b[((ihi)-1) + ((ihi - 2) - 1) * ldb];
+        h[(1 - 1) + (2 - 1) * 2] = b[((ihi - 1) - 1) + ((ihi - 1) - 1) * ldb];
+        h[(2 - 1) + (2 - 1) * 2] = b[((ihi)-1) + ((ihi - 1) - 1) * ldb];
+        h[(1 - 1) + (3 - 1) * 2] = b[((ihi - 1) - 1) + ((ihi)-1) * ldb];
+        h[(2 - 1) + (3 - 1) * 2] = b[((ihi)-1) + ((ihi)-1) * ldb];
         // Make H upper triangular
         Rlartg(h[0], h[(2 - 1)], c1, s1, temp);
         h[(2 - 1)] = zero;
@@ -97,7 +94,13 @@
         //
         // Normal operation, move bulge down
         //
-        h = b(k + 1, k + 2, k, k + 2);
+        // H = B( K+1:K+2, K:K+2 )  -- 2x3 submatrix
+        h[(1 - 1) + (1 - 1) * 2] = b[((k + 1) - 1) + ((k)-1) * ldb];
+        h[(2 - 1) + (1 - 1) * 2] = b[((k + 2) - 1) + ((k)-1) * ldb];
+        h[(1 - 1) + (2 - 1) * 2] = b[((k + 1) - 1) + ((k + 1) - 1) * ldb];
+        h[(2 - 1) + (2 - 1) * 2] = b[((k + 2) - 1) + ((k + 1) - 1) * ldb];
+        h[(1 - 1) + (3 - 1) * 2] = b[((k + 1) - 1) + ((k + 2) - 1) * ldb];
+        h[(2 - 1) + (3 - 1) * 2] = b[((k + 2) - 1) + ((k + 2) - 1) * ldb];
         //
         // Make H upper triangular
         //
