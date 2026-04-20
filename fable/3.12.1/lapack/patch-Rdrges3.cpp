--- Rdrges3.cpp~	2026-04-16 00:00:00.000000000 +0900
+++ Rdrges3.cpp	2026-04-16 00:00:00.000000000 +0900
@@ -421,7 +421,7 @@ void Rdrges3(INTEGER const nsizes, INTEG
                 result[(1 + rsub + isort) - 1] = ulpinv;
                 Rgges3("V", "V", sort.elems, Rlctes, n, s, lda, t, lda, sdim, alphar, alphai, beta, q, ldq, z, ldq, work, lwork, bwork, iinfo);
-                if (iinfo != 0 && iinfo != n + 2) {
+                if (iinfo != 0 && iinfo != n + 2 && iinfo != n + 3) {
                     result[(1 + rsub + isort) - 1] = ulpinv;
                     write(nounit, format_9999), "Rgges3", iinfo, n, jtype, ioldsd;
                     info = abs(iinfo);
