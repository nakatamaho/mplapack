diff --git a/mplapack/reference/Ctrsyl3.cpp b/mplapack/reference/Ctrsyl3.cpp
index 9346f9bb..b1ec31e9 100644
--- a/mplapack/reference/Ctrsyl3.cpp
+++ b/mplapack/reference/Ctrsyl3.cpp
@@ -38,24 +38,6 @@
 #include <memory>
 
 void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *c, INTEGER const ldc, REAL &scale, REAL *swork, INTEGER &ldswork, INTEGER &info) {
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
@@ -110,6 +92,7 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
     // Quick return if possible
     //
     const REAL one = 1.0;
+    const REAL two = 2.0;
     scale = one;
     if (m == 0 || n == 0) {
         return;
@@ -245,14 +228,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
                         // Mark the computation as pointless.
                         buf = zero;
                     } else {
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
@@ -277,14 +260,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -333,14 +316,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -416,14 +399,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -448,14 +431,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -503,14 +486,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -586,14 +569,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -618,14 +601,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -673,14 +656,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -756,14 +739,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -788,14 +771,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
@@ -844,14 +827,14 @@ void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER c
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
