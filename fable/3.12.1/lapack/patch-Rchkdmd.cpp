--- Rchkdmd.cpp
+++ Rchkdmd.cpp
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
@@ -426,7 +436,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                 lwork = 4 * m + 1;
                 work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
                 work = work_storage.get();
-                ac = a;
+                copy_matrix_block(ac, lda, 1, a, lda, 1, m, m);
                 // LAPACK CALL
                 Rgeev("N", "V", m, ac, m, reiga, ieiga, va, m, va, m, work, lwork, info);
                 tmp = zero;
@@ -452,22 +462,22 @@ void program_dmd_test(int argc, char const *argv[]) {
                     for (i = 1; i <= n / 2; i = i + 1) {
                         Rgemv("N", m, m, one, a, m, &f1[(i - 1) * ldf], 1, zero, &f1[((i + 1) - 1) * ldf], 1);
                     }
-                    x0(1, m, 1, n / 2) = f1(1, m, 1, n / 2);
-                    y0(1, m, 1, n / 2) = f1(1, m, 2, n / 2 + 1);
+                    copy_matrix_block(x0, ldx, 1, f1, ldf, 1, m, n / 2);
+                    copy_matrix_block(y0, ldy, 1, f1, ldf, 2, m, n / 2);
                     //
                     Rlarnv(2, iseed, m, &f1[0]);
                     for (i = 1; i <= n - n / 2; i = i + 1) {
                         Rgemv("N", m, m, one, a, m, &f1[(i - 1) * ldf], 1, zero, &f1[((i + 1) - 1) * ldf], 1);
                     }
-                    x0(1, m, n / 2 + 1, n) = f1(1, m, 1, n - n / 2);
-                    y0(1, m, n / 2 + 1, n) = f1(1, m, 2, n - n / 2 + 1);
+                    copy_matrix_block(x0, ldx, n / 2 + 1, f1, ldf, 1, m, n - n / 2);
+                    copy_matrix_block(y0, ldy, n / 2 + 1, f1, ldf, 2, m, n - n / 2);
                 } else {
                     Rlarnv(2, iseed, m, &f[0]);
                     for (i = 1; i <= n; i = i + 1) {
                         Rgemv("N", m, m, one, a, m, &f[(i - 1) * ldf], 1, zero, &f[((i + 1) - 1) * ldf], 1);
                     }
-                    x0(1, m, 1, n) = f(1, m, 1, n);
-                    y0(1, m, 1, n) = f(1, m, 2, n + 1);
+                    copy_matrix_block(x0, ldx, 1, f, ldf, 1, m, n);
+                    copy_matrix_block(y0, ldy, 1, f, ldf, 2, m, n);
                 }
                 //
                 xnorm = Rlange("F", m, n, x0, ldx, wdummy);
@@ -536,8 +546,8 @@ void program_dmd_test(int argc, char const *argv[]) {
                                         // Workspace query for the minimal (1) and for the optimal
                                         // (2) workspace lengths determined by workspace query.
                                         //
-                                        x(1, m, 1, n) = x0(1, m, 1, n);
-                                        y(1, m, 1, n) = y0(1, m, 1, n);
+                                        copy_matrix_block(x, ldx, 1, x0, ldx, 1, m, n);
+                                        copy_matrix_block(y, ldy, 1, y0, ldy, 1, m, n);
                                         //
                                         // Rgedmd: Workspace query and workspace allocation
                                         Rgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, x, ldx, y, ldy, nrnk, tol, k, reig, ieig, z, ldz, res, au, ldau, w, ldw, s, lds, wdummy, -1, idummy, -1, info);
@@ -553,7 +563,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                         Rgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, x, ldx, y, ldy, nrnk, tol, k, reig, ieig, z, ldz, res, au, ldau, w, ldw, s, lds, work, lwork, iwork, liwork, info);
                                         //
                                         for (INTEGER i_ = 1; i_ <= n; i_++) {
-                                            singvx[i_ - 1] = work[__SLICE__(1, n)];
+                                            singvx[i_ - 1] = work[i_ - 1];
                                         }
                                         //
                                         // ...... Rgedmd check point
@@ -706,7 +716,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                         // ======================================================================
                                         if (test_qrdmd && (k_traj == 1)) {
                                             rjobdata[2 - 1] = 1;
-                                            f1 = f;
+                                            copy_matrix_block(f1, ldf, 1, f, ldf, 1, m, n + 1);
                                             //
                                             // Rgedmdq test: Workspace query and workspace allocation
                                             Rgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, f1, ldf, x, ldx, y, ldy, nrnk, tol, kq, reigq, ieigq, z, ldz, res, au, ldau, w, ldw, s, lds, wdummy, -1, idummy, -1, info);
@@ -720,7 +730,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                             Rgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, f1, ldf, x, ldx, y, ldy, nrnk, tol, kq, reigq, ieigq, z, ldz, res, au, ldau, w, ldw, s, lds, work, lwork, iwork, liwork, info);
                                             //
                                             for (INTEGER i_ = 1; i_ <= kq; i_++) {
-                                                singvqx[i_ - 1] = work[__SLICE__(min(m, n + 1) + 1, min(m, n + 1) + kq)];
+                                                singvqx[i_ - 1] = work[min(m, n + 1) + i_ - 1];
                                             }
                                             //
                                             // ..... Rgedmdq check point
@@ -747,7 +757,7 @@ void program_dmd_test(int argc, char const *argv[]) {
                                                 // Check that the QR factors are computed and returned
                                                 // as requested. The residual ||F-Q*R||_F / ||F||_F
                                                 // is compared to M*N*EPS.
-                                                f2 = f;
+                                                copy_matrix_block(f2, ldf, 1, f, ldf, 1, m, n + 1);
                                                 Rgemm("N", "N", m, n + 1, min(m, n + 1), -one, f1, ldf, y, ldy, one, f2, ldf);
                                                 tmp_fqr = Rlange("F", m, n + 1, f2, ldf, work) / Rlange("F", m, n + 1, f, ldf, work);
                                                 if (tmp_fqr > tol2) {
--- Rchkdmd.cpp~	2026-03-25 13:28:22.828294940 +0900
+++ Rchkdmd.cpp	2026-03-25 13:29:26.049645010 +0900
@@ -298,11 +298,11 @@
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
@@ -748,7 +748,7 @@
                                                 nfail_svdiff++;
                                                 for (j = 1; j <= 3; j = j + 1) {
                                                     write(6, star), j, singvx[j - 1], singvqx[j - 1];
-                                                    read(6, star);
+                                                    read(5, star);
                                                 }
                                             }
                                             //
@@ -941,7 +941,6 @@
     //
     write(6, star);
     write(6, star), "Test completed.";
-    FEM_STOP(0);
 }
 
 int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dmd_test); }
