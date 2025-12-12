diff --git a/mplapack/reference/Rlaln2.cpp b/mplapack/reference/Rlaln2.cpp
index 7e2f56c3..e99de9c8 100644
--- a/mplapack/reference/Rlaln2.cpp
+++ b/mplapack/reference/Rlaln2.cpp
@@ -29,16 +29,23 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
-#define ci (equiv_0)
-#define cr (equiv_1)
-#define civ (equiv_0)
-#define crv (equiv_1)
-
-void Rlaln2(bool const ltrans, INTEGER const na, INTEGER const nw, REAL const smin, REAL const ca, REAL *a, INTEGER const lda, REAL const d1, REAL const d2, REAL *b, INTEGER const ldb, REAL const wr, REAL const wi, REAL *x, INTEGER const ldx, REAL &scale, REAL &xnorm, INTEGER &info) {
-    static bool zswap[] = {false, false, true, true};
-    static bool rswap[] = {false, true, false, true};
-    static INTEGER ipivot[] = {1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1};
+void Rlaln2(bool const ltrans, INTEGER const na, INTEGER const nw, REAL const &smin, REAL const &ca, REAL *a, INTEGER const lda, REAL const &d1, REAL const &d2, REAL *b, INTEGER const ldb, REAL const &wr, REAL const &wi, REAL *x, INTEGER const ldx, REAL &scale, REAL &xnorm, INTEGER &info) {
     INTEGER ipivot[4 * 4];
+    bool rswap[4];
+    bool zswap[4];
+    local_equivalences loc_equivalences;
+    {
+        using fem::mbr; // member
+        mbr<double> ci(dimension(2, 2));
+        mbr<double> civ(dimension(4));
+        mbr<double> cr(dimension(2, 2));
+        mbr<double> crv(dimension(4));
+        loc_equivalences.allocate(), equivalence(ci, civ).align<1>(arr_index(1, 1)).with<2>(arr_index(1)), equivalence(cr, crv).align<1>(arr_index(1, 1)).with<2>(arr_index(1));
+    }
+    arr_ref<double, 2> ci(loc_equivalences.bind<double>(), dimension(2, 2));
+    arr_ref<double> civ(loc_equivalences.bind<double>(), dimension(4));
+    arr_ref<double, 2> cr(loc_equivalences.bind<double>(), dimension(2, 2));
+    arr_ref<double> crv(loc_equivalences.bind<double>(), dimension(4));
     INTEGER ldcr = 2;
     INTEGER ldipivot = 4;
     INTEGER ldci = 2;
@@ -47,7 +54,6 @@ void Rlaln2(bool const ltrans, INTEGER const na, INTEGER const nw, REAL const sm
     //
     const REAL two = 2.0;
     REAL smlnum = two * Rlamch("Safe minimum");
-    // REAL smlnum = two * 2.2250738585072014E-308;
     const REAL one = 1.0;
     REAL bignum = one / smlnum;
     REAL smini = max(smin, smlnum);
@@ -95,8 +101,6 @@ void Rlaln2(bool const ltrans, INTEGER const na, INTEGER const nw, REAL const sm
     REAL bi1 = 0.0;
     REAL xi2 = 0.0;
     REAL xi1 = 0.0;
-    REAL equiv_0[4];
-    REAL equiv_1[4];
     if (na == 1) {
         //
         // 1 x 1  (i.e., scalar) system   C X = B
