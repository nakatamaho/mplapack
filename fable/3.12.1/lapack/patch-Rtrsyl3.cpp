--- a/mplapack/reference/Rtrsyl3.cpp
+++ b/mplapack/reference/Rtrsyl3.cpp
@@ -39,24 +39,6 @@
 
 void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL &scale, INTEGER *iwork, INTEGER const liwork, REAL *swork, INTEGER const ldswork, INTEGER &info) {
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
@@ -110,6 +92,7 @@
     // Quick return if possible
     //
     const REAL one = 1.0;
+    const REAL two = 2.0;
     scale = one;
     if (m == 0 || n == 0) {
         return;
@@ -315,14 +298,14 @@
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
@@ -347,14 +330,14 @@
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
@@ -403,14 +386,14 @@
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
@@ -486,14 +469,14 @@
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
@@ -518,14 +501,14 @@
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
@@ -573,14 +556,14 @@
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
@@ -657,14 +640,14 @@
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
@@ -688,14 +671,14 @@
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
@@ -743,14 +726,14 @@
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
@@ -826,14 +809,14 @@
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
@@ -858,14 +841,14 @@
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
@@ -914,14 +897,14 @@
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
