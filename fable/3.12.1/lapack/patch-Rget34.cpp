--- Rget34.cpp~	2026-04-15 00:00:00.000000000 +0900
+++ Rget34.cpp	2026-04-15 00:00:00.000000000 +0900
@@ -174,7 +174,13 @@ void Rget34(REAL &rmax, INTEGER &lmax, I
                                     if (t[(3 - 1) + (2 - 1) * ldt] != zero) {
                                         res += one / eps;
                                     }
-                                    if (t[(2 - 1)] != 0 && (t[0] != t[(2 - 1) + (2 - 1) * ldt] || sign(one, t[(2 - 1) * ldt]) == sign(one, t[(2 - 1)]))) {
+                                    bool diag12_mismatch = t[0] != t[(2 - 1) + (2 - 1) * ldt];
+#if defined ___MPLAPACK_BUILD_WITH_QD___
+                                    REAL diag12_diff = t[0] - t[(2 - 1) + (2 - 1) * ldt];
+                                    diag12_diff.renorm();
+                                    diag12_mismatch = !diag12_diff.is_zero();
+#endif
+                                    if (t[(2 - 1)] != 0 && (diag12_mismatch || sign(one, t[(2 - 1) * ldt]) == sign(one, t[(2 - 1)]))) {
                                         res += one / eps;
                                     }
                                 }
@@ -223,7 +229,13 @@ void Rget34(REAL &rmax, INTEGER &lmax, I
                                 Rhst01(3, 1, 3, t1, 4, t, 4, q, 4, work, lwork, result);
                                 res = result[1 - 1] + result[2 - 1];
                                 if (info == 0) {
-                                    if (t1[(3 - 1) + (3 - 1) * ldt1] != t[0]) {
+                                    bool diag11_mismatch = t1[(3 - 1) + (3 - 1) * ldt1] != t[0];
+#if defined ___MPLAPACK_BUILD_WITH_QD___
+                                    REAL diag11_diff = t1[(3 - 1) + (3 - 1) * ldt1] - t[0];
+                                    diag11_diff.renorm();
+                                    diag11_mismatch = !diag11_diff.is_zero();
+#endif
+                                    if (diag11_mismatch) {
                                         res += one / eps;
                                     }
                                     if (t[(2 - 1)] != zero) {
@@ -232,7 +244,13 @@ void Rget34(REAL &rmax, INTEGER &lmax, I
                                     if (t[(3 - 1)] != zero) {
                                         res += one / eps;
                                     }
-                                    if (t[(3 - 1) + (2 - 1) * ldt] != 0 && (t[(2 - 1) + (2 - 1) * ldt] != t[(3 - 1) + (3 - 1) * ldt] || sign(one, t[(2 - 1) + (3 - 1) * ldt]) == sign(one, t[(3 - 1) + (2 - 1) * ldt]))) {
+                                    bool diag23_mismatch = t[(2 - 1) + (2 - 1) * ldt] != t[(3 - 1) + (3 - 1) * ldt];
+#if defined ___MPLAPACK_BUILD_WITH_QD___
+                                    REAL diag23_diff = t[(2 - 1) + (2 - 1) * ldt] - t[(3 - 1) + (3 - 1) * ldt];
+                                    diag23_diff.renorm();
+                                    diag23_mismatch = !diag23_diff.is_zero();
+#endif
+                                    if (t[(3 - 1) + (2 - 1) * ldt] != 0 && (diag23_mismatch || sign(one, t[(2 - 1) + (3 - 1) * ldt]) == sign(one, t[(3 - 1) + (2 - 1) * ldt]))) {
                                         res += one / eps;
                                     }
                                 }
