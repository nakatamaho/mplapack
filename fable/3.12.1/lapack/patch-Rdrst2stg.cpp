https://claude.ai/chat/b5d92812-9d0f-45c3-827f-d03b184f7f89

diff --git a/mplapack/test/eig/common/Rdrvst2stg.cpp b/mplapack/test/eig/common/Rdrvst2stg.cpp
index 9401b8d..6abfe75 100644
--- Rdrvst2stg.cpp
+++ Rdrvst2stg.cpp
@@ -2150,7 +2150,7 @@ void Rdrvst2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *d
                         result[ntest - 1] = ulpinv;
                         result[(ntest + 1) - 1] = ulpinv;
                         result[(ntest + 2) - 1] = ulpinv;
-                        goto statement_700;
+                        goto statement_1720;
                     }
                 }
                 //
@@ -2171,13 +2171,13 @@ void Rdrvst2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *d
                         return;
                     } else {
                         result[ntest - 1] = ulpinv;
-                        goto statement_700;
+                        goto statement_1720;
                     }
                 }
                 //
                 if (m3 == 0 && n > 0) {
                     result[ntest - 1] = ulpinv;
-                    goto statement_700;
+                    goto statement_1720;
                 }
                 //
                 // Do test 78 (or +54)
@@ -2191,6 +2191,7 @@ void Rdrvst2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *d
                 }
                 result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
                 //
+            statement_1720:
                 Rlacpy(" ", n, n, v, ldu, a, lda);
                 //
             }
