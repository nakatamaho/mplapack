--- Cchkdmd.cpp
+++ Cchkdmd.cpp
@@ -44,6 +44,16 @@ using fem::common;
 #include <mplapack_eig.h>
 #include <memory>
 
+template <class T>
+static inline void copy_matrix_block(T *dst, INTEGER ldd, INTEGER dst_col, const T *src, INTEGER lds, INTEGER src_col, INTEGER rows, INTEGER cols) {
+    for (INTEGER j = 0; j < cols; j = j + 1) {
+        for (INTEGER i = 0; i < rows; i = i + 1) {
+            dst[i + (dst_col - 1 + j) * ldd] = src[i + (src_col - 1 + j) * lds];
+        }
+    }
+}
+
+
 // This is a test program for checking the implementations of
 // the implementations of the following subroutines
 //
@@ -377,7 +387,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                 zwork = zwork_storage.get();
                 work_storage = std::make_unique<REAL[]>(max((INTEGER)1, (2 * m)));
                 work = work_storage.get();
-                zac(1, m, 1, m) = za(1, m, 1, m);
+                copy_matrix_block(zac, lda, 1, za, lda, 1, m, m);
                 // LAPACK CALL
                 Cgeev("N", "N", m, zac, lda, zeigsa, zdum2x2, 2, zdum2x2, 2, zwork, lzwork, work, info);
                 //
@@ -395,23 +405,23 @@ void program_dmd_test(int argc, char const *argv[]) {
                     for (i = 1; i <= n / 2; i = i + 1) {
                         Cgemv("N", m, m, zone, za, lda, &zf[(i - 1) * ldf], 1, zzero, &zf[((i + 1) - 1) * ldf], 1);
                     }
-                    zx0(1, m, 1, n / 2) = zf(1, m, 1, n / 2);
-                    zy0(1, m, 1, n / 2) = zf(1, m, 2, n / 2 + 1);
+                    copy_matrix_block(zx0, ldx, 1, zf, ldf, 1, m, n / 2);
+                    copy_matrix_block(zy0, ldy, 1, zf, ldf, 2, m, n / 2);
                     //
                     Clarnv(2, iseed, m, &zf[0]);
                     for (i = 1; i <= n - n / 2; i = i + 1) {
                         Cgemv("N", m, m, zone, za, lda, &zf[(i - 1) * ldf], 1, zzero, &zf[((i + 1) - 1) * ldf], 1);
                     }
-                    zx0(1, m, n / 2 + 1, n) = zf(1, m, 1, n - n / 2);
-                    zy0(1, m, n / 2 + 1, n) = zf(1, m, 2, n - n / 2 + 1);
+                    copy_matrix_block(zx0, ldx, n / 2 + 1, zf, ldf, 1, m, n - n / 2);
+                    copy_matrix_block(zy0, ldy, n / 2 + 1, zf, ldf, 2, m, n - n / 2);
                 } else {
                     Clarnv(2, iseed, m, &zf[0]);
                     for (i = 1; i <= n; i = i + 1) {
                         Cgemv("N", m, m, zone, za, m, &zf[(i - 1) * ldf], 1, zzero, &zf[((i + 1) - 1) * ldf], 1);
                     }
-                    zf0(1, m, 1, n + 1) = zf(1, m, 1, n + 1);
-                    zx0(1, m, 1, n) = zf0(1, m, 1, n);
-                    zy0(1, m, 1, n) = zf0(1, m, 2, n + 1);
+                    copy_matrix_block(zf0, ldf, 1, zf, ldf, 1, m, n + 1);
+                    copy_matrix_block(zx0, ldx, 1, zf0, ldf, 1, m, n);
+                    copy_matrix_block(zy0, ldy, 1, zf0, ldf, 2, m, n);
                 }
                 //
                 // ........................................................................
@@ -471,8 +481,8 @@ void program_dmd_test(int argc, char const *argv[]) {
                                         // Cgedmd is always tested and its results are also used for
                                         // comparisons with Cgedmdq.
                                         //
-                                        zx(1, m, 1, n) = zx0(1, m, 1, n);
-                                        zy(1, m, 1, n) = zy0(1, m, 1, n);
+                                        copy_matrix_block(zx, ldx, 1, zx0, ldx, 1, m, n);
+                                        copy_matrix_block(zy, ldy, 1, zy0, ldy, 1, m, n);
                                         //
                                         Cgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, zx, ldx, zy, ldy, nrnk, tol, k, zeigs, zz, ldz, res, zau, ldau, zw, ldw, zs, lds, zdummy, -1, wdummy, -1, idummy, -1, info);
                                         if ((info == 2) || (info == 3) || (info < 0)) {
@@ -505,7 +515,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                         }
                                         //
                                         for (INTEGER i_ = 1; i_ <= n; i_++) {
-                                            singvx[i_ - 1] = work[__SLICE__(1, n)];
+                                            singvx[i_ - 1] = work[i_ - 1];
                                         }
                                         //
                                         // ...... Cgedmd check point
@@ -611,7 +621,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                         //
                                         if (test_qrdmd && (k_traj == 1)) {
                                             //
-                                            zf(1, m, 1, n + 1) = zf0(1, m, 1, n + 1);
+                                            copy_matrix_block(zf, ldf, 1, zf0, ldf, 1, m, n + 1);
                                             //
                                             Cgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, zf, ldf, zx, ldx, zy, ldy, nrnk, tol, k, zeigs, zz, ldz, res, zau, ldau, zw, ldw, zs, lds, zdummy, -1, wdummy, -1, idummy, -1, info);
                                             //
@@ -635,7 +645,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                                 FEM_STOP(0);
                                             }
                                             for (INTEGER i_ = 1; i_ <= n; i_++) {
-                                                singvqx[i_ - 1] = work[__SLICE__(1, n)];
+                                                singvqx[i_ - 1] = work[i_ - 1];
                                             }
                                             //
                                             // ..... Cgedmdq check point
@@ -663,7 +673,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                                 // Check that the QR factors are computed and returned
                                                 // as requested. The residual ||F-Q*R||_F / ||F||_F
                                                 // is compared to M*N*EPS.
-                                                zf1(1, m, 1, n + 1) = zf0(1, m, 1, n + 1);
+                                                copy_matrix_block(zf1, ldf, 1, zf0, ldf, 1, m, n + 1);
                                                 Cgemm("N", "N", m, n + 1, min(m, n + 1), -zone, zf, ldf, zy, ldy, zone, zf1, ldf);
                                                 tmp_fqr = Clange("F", m, n + 1, zf1, ldf, work) / Clange("F", m, n + 1, zf0, ldf, work);
                                                 if (tmp_fqr > tol2) {
--- Cchkdmd.cpp~	2026-03-25 13:28:22.787294063 +0900
+++ Cchkdmd.cpp	2026-03-25 13:29:08.516270861 +0900
@@ -268,11 +268,11 @@
         //
         // Set the dimensions of the problem ...
         write(6, star), "M = ";
-        read(6, star), m;
+        read(5, star), m;
         write(6, star), m;
         // ... and the number of snapshots.
         write(6, star), "N = ";
-        read(6, star), n;
+        read(5, star), n;
         write(6, star), n;
         //
         // ... Test the dimensions
@@ -662,7 +662,7 @@
                                                     nfail_svdiff++;
                                                     for (j = 1; j <= 3; j = j + 1) {
                                                         write(6, star), j, singvx[j - 1], singvqx[j - 1];
-                                                        read(6, star);
+                                                        read(5, star);
                                                     }
                                                     //
                                                 }
@@ -832,7 +832,6 @@
     //
     write(6, star);
     write(6, star), "Test completed.";
-    FEM_STOP(0);
 }
 
 int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dmd_test); }
