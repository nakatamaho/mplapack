diff --git a/mplapack/reference/Rtrsyl3.cpp b/mplapack/reference/Rtrsyl3.cpp
index e8d9ca25..a5df9f29 100644
--- a/mplapack/reference/Rtrsyl3.cpp
+++ b/mplapack/reference/Rtrsyl3.cpp
@@ -38,24 +38,6 @@
 #include <memory>
 
 void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL &scale, INTEGER *iwork, INTEGER const liwork, REAL *swork, INTEGER &ldswork, INTEGER &info) {
-    //
-    // .. Scalar Arguments ..
-    // ..
-    // .. Array Arguments ..
-    // ..
-    // .. Parameters ..
-    // ..
-    // .. Local Scalars ..
-    // ..
-    // .. Local Arrays ..
-    // ..
-    // .. External Functions ..
-    // ..
-    // .. External Subroutines ..
-    // ..
-    // .. Intrinsic Functions ..
-    // ..
-    // .. Executable Statements ..
     //
     // Decode and Test input parameters
     //
@@ -111,6 +93,7 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
     // Quick return if possible
     //
     const REAL one = 1.0;
+    const REAL two = 2.0;
     scale = one;
     if (m == 0 || n == 0) {
         return;
@@ -316,14 +299,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                         buf = zero;
                     } else {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                     }
                     for (jj = 1; jj <= nbb; jj = jj + 1) {
                         for (ll = 1; ll <= nba; ll = ll + 1) {
                             // Bound by BIGNUM to not introduce Inf. The value
                             // is irrelevant; corresponding entries of the
                             // solution will be flushed in consistency scaling.
-                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                         }
                     }
                 }
@@ -348,14 +331,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(anrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -404,14 +387,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(bnrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -487,14 +470,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                         buf = zero;
                     } else {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                     }
                     for (jj = 1; jj <= nbb; jj = jj + 1) {
                         for (ll = 1; ll <= nba; ll = ll + 1) {
                             // Bound by BIGNUM to not introduce Inf. The value
                             // is irrelevant; corresponding entries of the
                             // solution will be flushed in consistency scaling.
-                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                         }
                     }
                 }
@@ -519,14 +502,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(anrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -574,14 +557,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(bnrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -658,14 +641,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                         buf = zero;
                     } else {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                     }
                     for (jj = 1; jj <= nbb; jj = jj + 1) {
                         for (ll = 1; ll <= nba; ll = ll + 1) {
                             // Bound by BIGNUM to not introduce Inf. The value
                             // is irrelevant; corresponding entries of the
                             // solution will be flushed in consistency scaling.
-                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                         }
                     }
                 }
@@ -689,14 +672,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(anrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -744,14 +727,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(bnrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -827,14 +810,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                         buf = zero;
                     } else {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                     }
                     for (jj = 1; jj <= nbb; jj = jj + 1) {
                         for (ll = 1; ll <= nba; ll = ll + 1) {
                             // Bound by BIGNUM to not introduce Inf. The value
                             // is irrelevant; corresponding entries of the
                             // solution will be flushed in consistency scaling.
-                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                         }
                     }
                 }
@@ -859,14 +842,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(anrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
@@ -915,14 +898,14 @@ void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                     scaloc = Rlarmm(bnrm, xnrm, cnrm);
                     if (scaloc * scamin == zero) {
                         // Use second scaling factor to prevent flushing to zero.
-                        buf = buf * pow(2.0, Mexponent(scaloc));
+                        buf = buf * pow(two, Mexponent(scaloc));
                         for (jj = 1; jj <= nbb; jj = jj + 1) {
                             for (ll = 1; ll <= nba; ll = ll + 1) {
-                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(2.0, Mexponent(scaloc)));
+                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                             }
                         }
-                        scamin = scamin / pow(2.0, Mexponent(scaloc));
-                        scaloc = scaloc / pow(2.0, Mexponent(scaloc));
+                        scamin = scamin / pow(two, Mexponent(scaloc));
+                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                     }
                     cnrm = cnrm * scaloc;
                     xnrm = xnrm * scaloc;
