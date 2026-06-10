--- Rdrvsg.cpp
+++ Rdrvsg.cpp
@@ -329,8 +329,8 @@ void Rdrvsg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotyp
                 il = 1;
                 iu = n;
             } else {
-                il = 1 + (n - 1) * Rlarnd(1, iseed2);
-                iu = 1 + (n - 1) * Rlarnd(1, iseed2);
+                il = (INTEGER)1 + castINTEGER((n - 1) * Rlarnd(1, iseed2));
+                iu = (INTEGER)1 + castINTEGER((n - 1) * Rlarnd(1, iseed2));
                 if (il > iu) {
                     itemp = il;
                     il = iu;
