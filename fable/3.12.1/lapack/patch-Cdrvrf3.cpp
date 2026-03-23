diff --git a/mplapack/test/lin/common/Cdrvrf3.cpp b/mplapack/test/lin/common/Cdrvrf3.cpp
index 3ae5864a..e23867b1 100644
--- Cdrvrf3.cpp
+++ Cdrvrf3.cpp
@@ -96,6 +96,7 @@ void Cdrvrf3(INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL const thr
     INTEGER j = 0;
     const INTEGER ntests = 1;
     REAL result[ntests];
+    const REAL two = 2.0;
     for (iim = 1; iim <= nn; iim = iim + 1) {
         //
         m = nval[iim - 1];
@@ -178,6 +179,19 @@ void Cdrvrf3(INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL const thr
                                         //
                                         srnamt = "Cgeqrf";
                                         Cgeqrf(na, na, a, lda, tau, z_work_zgeqrf, lda, info);
+                                        //
+                                        // Forcing main diagonal of test matrix to
+                                        // be unit makes it ill-conditioned for
+                                        // some test cases
+                                        //
+                                        if (Mlsame(diag.elems, "U")) {
+                                            for (j = 1; j <= na; j = j + 1) {
+                                                for (i = 1; i <= j; i = i + 1) {
+                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] / (two * a[(j - 1) + (j - 1) * lda]);
+                                                }
+                                            }
+                                        }
+                                        //
                                     } else {
                                         //
                                         // The case IUPLO.EQ.2 is when SIDE.EQ.'L'
@@ -185,6 +199,19 @@ void Cdrvrf3(INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL const thr
                                         //
                                         srnamt = "Cgelqf";
                                         Cgelqf(na, na, a, lda, tau, z_work_zgeqrf, lda, info);
+                                        //
+                                        // Forcing main diagonal of test matrix to
+                                        // be unit makes it ill-conditioned for
+                                        // some test cases
+                                        //
+                                        if (Mlsame(diag.elems, "U")) {
+                                            for (i = 1; i <= na; i = i + 1) {
+                                                for (j = 1; j <= i; j = j + 1) {
+                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] / (two * a[(i - 1) + (i - 1) * lda]);
+                                                }
+                                            }
+                                        }
+                                        //
                                     }
                                     //
                                     // After the QR factorization, the diagonal
