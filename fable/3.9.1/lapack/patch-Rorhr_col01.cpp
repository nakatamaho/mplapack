--- Rorhr_col01.cpp
+++ Rorhr_col01.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Rorhr_col01(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -52,7 +54,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     //
     // Dynamically allocate local arrays
     //
@@ -81,8 +83,6 @@
     //
     INTEGER nrb = max((INTEGER)1, iceil(castREAL(m - n) / castREAL(mb1 - n)));
     //
-    // FABLE: ALLOCATE removed (RAII in C++)
-    //
     // Begin determine LWORK for the array WORK and allocate memory.
     //
     // Rlatsqr requires NB1 to be bounded by N.
@@ -298,11 +298,6 @@
         result[6 - 1] = zero;
     }
     //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t1,t2,diag,c,d,cf,d"
-                        "f)");
-    //
     // End of Rorhr_col01
     //
 }
