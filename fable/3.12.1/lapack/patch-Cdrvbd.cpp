--- Cdrvbd.cpp~	2026-04-16 00:00:00.000000000 +0900
+++ Cdrvbd.cpp	2026-04-16 00:00:00.000000000 +0900
@@ -598,6 +598,11 @@
                     Clacpy("F", m, n, asav, lda, usav, lda);
                     srnamt = "Cgesvj";
                     Cgesvj("G", "U", "V", m, n, usav, lda, ssav, 0, a, ldvt, work, lwork, rwork, lrwork, iinfo);
+                    if (rwork[1 - 1] != one) {
+                        for (i = 1; i <= mnmin; i = i + 1) {
+                            ssav[i - 1] = rwork[1 - 1] * ssav[i - 1];
+                        }
+                    }
                     //
                     // Cgesvj returns V not VH
                     //
