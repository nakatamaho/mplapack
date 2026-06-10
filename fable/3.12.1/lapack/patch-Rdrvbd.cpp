--- Rdrvbd.cpp~	2026-04-16 00:00:00.000000000 +0900
+++ Rdrvbd.cpp	2026-04-16 00:00:00.000000000 +0900
@@ -572,6 +572,11 @@
                     Rlacpy("F", m, n, asav, lda, usav, lda);
                     srnamt = "Rgesvj";
                     Rgesvj("G", "U", "V", m, n, usav, lda, ssav, 0, a, ldvt, work, lwork, info);
+                    if (work[1 - 1] != one) {
+                        for (i = 1; i <= mnmin; i = i + 1) {
+                            ssav[i - 1] = work[1 - 1] * ssav[i - 1];
+                        }
+                    }
                     //
                     // Rgesvj returns V not VT
                     //
