--- a/mplapack/reference/Ctrsyl3.cpp
+++ b/mplapack/reference/Ctrsyl3.cpp
@@ -39,24 +39,6 @@
 
 void Ctrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *c, INTEGER const ldc, REAL &scale, REAL *swork, INTEGER const ldswork, INTEGER &info) {
     //
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
-    //
     // Decode and Test input parameters
     //
     bool notrna = Mlsame(trana, "N");
@@ -109,6 +91,7 @@
     // Quick return if possible
     //
     const REAL one = 1.0;
+    const REAL two = 2.0;
     scale = one;
     if (m == 0 || n == 0) {
         return;
@@ -244,14 +227,14 @@
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
@@ -276,14 +259,14 @@
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
@@ -332,14 +315,14 @@
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
@@ -415,14 +398,14 @@
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
@@ -447,14 +430,14 @@
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
@@ -502,14 +485,14 @@
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
@@ -585,14 +568,14 @@
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
@@ -617,14 +600,14 @@
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
@@ -672,14 +655,14 @@
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
@@ -755,14 +738,14 @@
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
@@ -787,14 +770,14 @@
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
@@ -843,14 +826,14 @@
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
