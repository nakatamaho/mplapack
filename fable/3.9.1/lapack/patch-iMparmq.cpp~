diff --git a/mplapack/reference/Rbdsvdx.cpp b/mplapack/reference/Rbdsvdx.cpp
index 43876710..7ce454c9 100644
--- a/mplapack/reference/Rbdsvdx.cpp
+++ b/mplapack/reference/Rbdsvdx.cpp
@@ -29,7 +29,7 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER const n, REAL *d, REAL *e, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, INTEGER &ns, REAL *s, REAL *z, INTEGER const ldz, REAL *work, INTEGER *iwork, INTEGER &info) {
+void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER const n, REAL *d, REAL *e, REAL const &vl, REAL const &vu, INTEGER const il, INTEGER const iu, INTEGER &ns, REAL *s, REAL *z, INTEGER const ldz, REAL *work, INTEGER *iwork, INTEGER &info) {
     //
     // Test the input parameters.
     //
@@ -100,21 +100,20 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
         return;
     }
     //
-    REAL two = 2.0;
-    REAL abstol = two * Rlamch("Safe Minimum");
+    REAL abstol = 2 * Rlamch("Safe Minimum");
     REAL ulp = Rlamch("Precision");
     REAL eps = Rlamch("Epsilon");
-    REAL sqrt2 = sqrt(two);
+    REAL sqrt2 = sqrt(2.0);
     REAL ortol = sqrt(ulp);
     //
     // Criterion for splitting is taken from Rbdsqr when singular
     // values are computed to relative accuracy TOL. (See J. Demmel and
     // W. Kahan, Accurate singular values of bidiagonal matrices, SIAM
-    // J. Sci. and Stat. Comput., 11:873912, 1990.)
+    // J. Sci. and Stat. Comput., 11:873–912, 1990.)
     //
     const REAL ten = 10.0;
     const REAL hndrd = 100.0;
-    const REAL meigth = -0.1250e0;
+    const REAL meigth = -0.125;
     REAL tol = max(ten, min(hndrd, pow(eps, meigth))) * eps;
     //
     // Compute approximate maximum, minimum singular values.
@@ -122,7 +121,7 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
     INTEGER i = iRamax(n, d, 1);
     REAL smax = abs(d[i - 1]);
     i = iRamax(n - 1, e, 1);
-    smax = max(smax, REAL(abs(e[i - 1])));
+    smax = max(smax, abs(e[i - 1]));
     //
     // Compute threshold for neglecting D's and E's.
     //
@@ -178,6 +177,7 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
         //
         // All singular values will be found. We aim at -s (see
         // leading comments) with RNGVX = 'I'. IL and IU are set
+        // later (as ILTGK and IUTGK) according to the dimension
         // of the active submatrix.
         //
         rngvx = 'I';
@@ -193,8 +193,9 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
         rngvx = 'V';
         vltgk = -vu;
         vutgk = -vl;
-        for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
-            work[l - 1] = zero;
+        for (INTEGER i = idtgk; i <= idtgk + 2 * n - 1; i++) {
+            work[i - 1] = zero;
+        }
         Rcopy(n, d, 1, &work[ietgk - 1], 2);
         Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
         Rstevx("N", "V", n * 2, &work[idtgk - 1], &work[ietgk - 1], vltgk, vutgk, iltgk, iltgk, abstol, ns, s, z, ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
@@ -217,14 +218,16 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
         iltgk = il;
         iutgk = iu;
         rngvx = 'V';
-        for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
-            work[l - 1] = zero;
+        for (INTEGER i = idtgk; i <= idtgk + 2 * n - 1; i++) {
+            work[i - 1] = zero;
+        }
         Rcopy(n, d, 1, &work[ietgk - 1], 2);
         Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
         Rstevx("N", "I", n * 2, &work[idtgk - 1], &work[ietgk - 1], vltgk, vltgk, iltgk, iltgk, abstol, ns, s, z, ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
         vltgk = s[0] - fudge * smax * ulp * n;
-        for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
-            work[l - 1] = zero;
+        for (INTEGER i = idtgk; i <= idtgk + 2 * n - 1; i++) {
+            work[i - 1] = zero;
+        }
         Rcopy(n, d, 1, &work[ietgk - 1], 2);
         Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
         Rstevx("N", "I", n * 2, &work[idtgk - 1], &work[ietgk - 1], vutgk, vutgk, iutgk, iutgk, abstol, ns, s, z, ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
@@ -264,11 +267,13 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
     //
     // Form the tridiagonal TGK matrix.
     //
-    for (int i = 1; i <= n; i++)
+    for (INTEGER i = 1; i <= n; i++) {
         s[i - 1] = zero;
+    }
     work[(ietgk + 2 * n - 1) - 1] = zero;
-    for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
-        work[l - 1] = zero;
+    for (INTEGER i = idtgk; i <= idtgk + 2 * n - 1; i++) {
+        work[i - 1] = zero;
+    }
     Rcopy(n, d, 1, &work[ietgk - 1], 2);
     Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
     //
@@ -382,7 +387,7 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
                         // Exit with the error code from Rstevx.
                         return;
                     }
-                    emin = abs(Mmaxval(s, isbeg, isbeg + nsl - 1, 1));
+                    emin = abs(maxval(s[__SLICE__(isbeg, isbeg + nsl - 1)]));
                     //
                     if (nsl > 0 && wantz) {
                         //
@@ -401,10 +406,8 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
                             // eigenvectors corresponding to the two smallest
                             // eigenvalues.
                             //
-                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
-                                z[(l - 1) + ((icolz + nsl - 2) - 1) * ldz] = z[(l - 1) + ((icolz + nsl - 2) - 1) * ldz] + z[(l - 1) + ((icolz + nsl - 1) - 1) * ldz];
-                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
-                                z[(l - 1) + ((icolz + nsl - 1) - 1) * ldz] = 0.0;
+                            z[__SLICE2D__(irowz, irowz + ntgk - 1, icolz + nsl - 2, ldz)] += z[__SLICE2D__(irowz, irowz + ntgk - 1, icolz + nsl - 1, ldz)];
+                            z[__SLICE2D__(irowz, irowz + ntgk - 1, icolz + nsl - 1, ldz)] = zero;
                             // IF( IUTGK*2.GT.NTGK ) THEN
                             // Eigenvalue equal to zero or very small.
                             // NSL = NSL - 1
@@ -451,11 +454,10 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
                             // active submatrix is reached).
                             //
                             split = true;
-
-                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
-                                z[l + ((n + 1) - 1) * ldz] = z[l + ((ns + nsl) - 1) * ldz];
-                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
-                                z[l + ((ns + nsl) - 1) * ldz] = 0.0;
+                            z[__SLICE2D__(irowz, irowz + ntgk - 1, n + 1, ldz)] = z[__SLICE2D__(irowz, irowz + ntgk - 1, ns + nsl, ldz)];
+                            for (INTEGER i = irowz; i <= irowz + ntgk - 1; i++) {
+                                z[(i - 1) + ((ns + nsl) - 1) * ldz] = zero;
+                            }
                         }
                         // ** WANTZ **!
                     }
@@ -483,8 +485,9 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
                     // ** NTGK.GT.0 **!
                 }
                 if (irowz < n * 2 && wantz) {
-                    for (int l = 1; l <= irowz - 1; l++)
-                        z[(l - 1) + (icolz - 1) * ldz] = zero;
+                    for (INTEGER i = 1; i <= irowz - 1; i++) {
+                        z[(i - 1) + (icolz - 1) * ldz] = zero;
+                    }
                 }
                 // ** IDPTR loop **!
             }
@@ -493,10 +496,10 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
                 // Bring back eigenvector corresponding
                 // to eigenvalue equal to zero.
                 //
-                for (int l = idbeg; l <= idend - ntgk + 1; l++)
-                    z[(l - 1) + (isbeg - 1) * ldz] = z[(l - 1) + (isbeg - 1) * ldz] + z[(l - 1) + ((n + 1) - 1) * ldz];
-                for (int l = idbeg; l <= idend - ntgk + 1; l++)
-                    z[(l - 1) + ((n + 1) - 1) * ldz] = 0.0;
+                z[__SLICE2D__(idbeg, idend - ntgk + 1, isbeg - 1, ldz)] += z[__SLICE2D__(idbeg, idend - ntgk + 1, n + 1, ldz)];
+                for (INTEGER i = idbeg; i <= idend - ntgk + 1; i++) {
+                    z[(i - 1) + ((n + 1) - 1) * ldz] = 0.0;
+                }
             }
             irowv = irowv - 1;
             irowu++;
@@ -535,12 +538,11 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
     if (indsv) {
         k = iu - il + 1;
         if (k < ns) {
-            for (int l = k + 1; l <= ns; l++)
-                s[l - 1] = zero;
+            for (INTEGER i = k + 1; i <= ns; i++) {
+                s[i - 1] = zero;
+            }
             if (wantz) {
-                for (int l = 1; l <= n * 2; l++)
-                    for (int m = k + 1; m <= ns; m++)
-                        z[(l - 1) + (m - 1) * ldz] = zero;
+                z(1, n * 2, k + 1, ns) = zero;
             }
             ns = k;
         }
@@ -553,10 +555,10 @@ void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER cons
         for (i = 1; i <= ns; i = i + 1) {
             Rcopy(n * 2, &z[(i - 1) * ldz], 1, work, 1);
             if (lower) {
-                Rcopy(n, &work[2 - 1], 2, &z[((n + 1) - 1) + (i - 1) * ldz], 1);
+                Rcopy(n, &work[1], 2, &z[((n + 1) - 1) + (i - 1) * ldz], 1);
                 Rcopy(n, &work[0], 2, &z[(i - 1) * ldz], 1);
             } else {
-                Rcopy(n, &work[2 - 1], 2, &z[(i - 1) * ldz], 1);
+                Rcopy(n, &work[1], 2, &z[(i - 1) * ldz], 1);
                 Rcopy(n, &work[0], 2, &z[((n + 1) - 1) + (i - 1) * ldz], 1);
             }
         }
diff --git a/mplapack/reference/iMparmq.cpp b/mplapack/reference/iMparmq.cpp
index 63917d6e..981f20b2 100644
--- a/mplapack/reference/iMparmq.cpp
+++ b/mplapack/reference/iMparmq.cpp
@@ -28,9 +28,6 @@
 
 #include <mpblas.h>
 #include <mplapack.h>
-#include <string.h>
-
-#define subnamlen 32
 
 INTEGER
 iMparmq(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER const /* n */, INTEGER const ilo, INTEGER const ihi, INTEGER const /* lwork */) {
@@ -40,7 +37,6 @@ iMparmq(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER
     const INTEGER iacc22 = 16;
     INTEGER nh = 0;
     INTEGER ns = 0;
-    INTEGER name_len;
     const REAL two = 2.0;
     if ((ispec == ishfts) || (ispec == inwin) || (ispec == iacc22)) {
         //
@@ -74,7 +70,7 @@ iMparmq(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER
     const INTEGER inibl = 14;
     const INTEGER nibble = 14;
     const INTEGER knwswp = 500;
-    char subnam[subnamlen];
+    char subnam[6];
     INTEGER ic = 0;
     INTEGER iz = 0;
     INTEGER i = 0;
@@ -124,19 +120,47 @@ iMparmq(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER
         // Convert NAME to upper case if the first character is lower case.
         //
         return_value = 0;
-        strncpy(subnam, name, subnamlen - 1);
-        ic = *subnam;
-        iz = 'Z';
+        subnam = name;
+        ic = fem::ichar(subnam[0]);
+        iz = fem::ichar("Z");
         if (iz == 90 || iz == 122) {
             //
             // ASCII character set
             //
             if (ic >= 97 && ic <= 122) {
-                *subnam = (char)(ic - 32);
-                for (i = 2; i <= 6; i++) {
-                    ic = subnam[i - 1];
+                subnam[0] = fem::fchar(ic - 32);
+                for (i = 2; i <= 6; i = i + 1) {
+                    ic = fem::ichar(subnam(i, i));
                     if (ic >= 97 && ic <= 122) {
-                        subnam[i - 1] = (char)(ic - 32);
+                        subnam(i, i) = fem::fchar(ic - 32);
+                    }
+                }
+            }
+            //
+        } else if (iz == 233 || iz == 169) {
+            //
+            // EBCDIC character set
+            //
+            if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
+                subnam[0] = fem::fchar(ic + 64);
+                for (i = 2; i <= 6; i = i + 1) {
+                    ic = fem::ichar(subnam(i, i));
+                    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
+                        subnam(i, i) = fem::fchar(ic + 64);
+                    }
+                }
+            }
+            //
+        } else if (iz == 218 || iz == 250) {
+            //
+            // Prime machines:  ASCII+128
+            //
+            if (ic >= 225 && ic <= 250) {
+                subnam[0] = fem::fchar(ic - 32);
+                for (i = 2; i <= 6; i = i + 1) {
+                    ic = fem::ichar(subnam(i, i));
+                    if (ic >= 225 && ic <= 250) {
+                        subnam(i, i) = fem::fchar(ic - 32);
                     }
                 }
             }
