--- Cunhr_col02.cpp
+++ Cunhr_col02.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Cunhr_col02(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, REAL *result) {
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
@@ -156,7 +158,7 @@
     // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
     //
     Claset("Full", m, m, czero, cone, r, m);
-    Cherk("U", "C", m, m, -cone, q, m, cone, r, m);
+    Cherk("U", "C", m, m, (-cone).real(), q, m, cone.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
@@ -256,11 +258,6 @@
         result[6 - 1] = zero;
     }
     //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t1,t2,diag,c,d,cf,d"
-                        "f)");
-    //
     // End of Cunhr_col02
     //
 }
