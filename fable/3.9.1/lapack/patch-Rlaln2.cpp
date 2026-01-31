--- a/mplapack/reference/Rlaln2.cpp
+++ b/mplapack/reference/Rlaln2.cpp
@@ -36,23 +36,15 @@
 #include <mpblas.h>
 #include <mplapack.h>
 
+#define ci (equiv_0)
+#define cr (equiv_1)
+#define civ (equiv_0)
+#define crv (equiv_1)
+
 void Rlaln2(bool const ltrans, INTEGER const na, INTEGER const nw, REAL const smin, REAL const ca, REAL *a, INTEGER const lda, REAL const d1, REAL const d2, REAL *b, INTEGER const ldb, REAL const wr, REAL const wi, REAL *x, INTEGER const ldx, REAL &scale, REAL &xnorm, INTEGER &info) {
-    static bool zswap[4] = {false, false, true, true};
-    static bool rswap[4] = {false, true, false, true};
+    static bool zswap[] = {false, false, true, true};
+    static bool rswap[] = {false, true, false, true};
     static INTEGER ipivot[] = {1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1};
-    local_equivalences loc_equivalences;
-    {
-        using fem::mbr; // member
-        mbr<double> ci(dimension(2, 2));
-        mbr<double> civ(dimension(4));
-        mbr<double> cr(dimension(2, 2));
-        mbr<double> crv(dimension(4));
-        loc_equivalences.allocate(), equivalence(ci, civ).align<1>(arr_index(1, 1)).with<2>(arr_index(1)), equivalence(cr, crv).align<1>(arr_index(1, 1)).with<2>(arr_index(1));
-    }
-    arr_ref<double, 2> ci(loc_equivalences.bind<double>(), dimension(2, 2));
-    arr_ref<double> civ(loc_equivalences.bind<double>(), dimension(4));
-    arr_ref<double, 2> cr(loc_equivalences.bind<double>(), dimension(2, 2));
-    arr_ref<double> crv(loc_equivalences.bind<double>(), dimension(4));
     INTEGER ldcr = 2;
     INTEGER ldipivot = 4;
     INTEGER ldci = 2;
@@ -61,6 +53,7 @@
     //
     const REAL two = 2.0;
     REAL smlnum = two * Rlamch("Safe minimum");
+    // REAL smlnum = two * 2.2250738585072014E-308;
     const REAL one = 1.0;
     REAL bignum = one / smlnum;
     REAL smini = max(smin, smlnum);
@@ -108,6 +101,8 @@
     REAL bi1 = 0.0;
     REAL xi2 = 0.0;
     REAL xi1 = 0.0;
+    REAL equiv_0[4];
+    REAL equiv_1[4];
     if (na == 1) {
         //
         // 1 x 1  (i.e., scalar) system   C X = B
@@ -258,7 +253,7 @@
             //
             xr2 = (br2 * scale) / ur22;
             xr1 = (scale * br1) * ur11r - xr2 * (ur11r * ur12);
-            if (Cswap[icmax - 1]) {
+            if (zswap[icmax - 1]) {
                 x[0] = xr2;
                 x[(2 - 1)] = xr1;
             } else {
@@ -395,7 +390,7 @@
             Rladiv(br2, bi2, ur22, ui22, xr2, xi2);
             xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
             xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
-            if (Cswap[icmax - 1]) {
+            if (zswap[icmax - 1]) {
                 x[0] = xr2;
                 x[(2 - 1)] = xr1;
                 x[(2 - 1) * ldx] = xi2;
