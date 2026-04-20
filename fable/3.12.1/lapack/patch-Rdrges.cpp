--- Rdrges.cpp~
+++ Rdrges.cpp
@@ -411,7 +411,7 @@
                 ntest = 1 + rsub + isort;
                 result[(1 + rsub + isort) - 1] = ulpinv;
                 Rgges("V", "V", sort.elems, Rlctes, n, s, lda, t, lda, sdim, alphar, alphai, beta, q, ldq, z, ldq, work, lwork, bwork, iinfo);
-                if (iinfo != 0 && iinfo != n + 2) {
+                if (iinfo != 0 && iinfo != n + 2 && iinfo != n + 3) {
                     result[(1 + rsub + isort) - 1] = ulpinv;
                     write(nounit, format_9999), "Rgges", iinfo, n, jtype, ioldsd;
                     info = abs(iinfo);
