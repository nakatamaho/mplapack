--- Cunhr_col01.cpp
+++ Cunhr_col01.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Cunhr_col01(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -52,11 +54,7 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
-    //
-    // Dynamically allocate local arrays
-    //
-    // FABLE: ALLOCATE removed (RAII in C++)
+    INTEGER l = max(m, n, (INTEGER)1);
     //
     // Put random numbers into A and copy to AF
     //
@@ -199,7 +197,7 @@
     // Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)
     //
     Claset("Full", m, m, czero, cone, r, m);
-    Cherk("U", "C", m, m, -cone, q, m, cone, r, m);
+    Cherk("U", "C", m, m, (-cone).real(), q, m, cone.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
@ -297,6 +297,11 @@ void Cunhr_col01(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER co
         result[6 - 1] = zero;
				     }
     //
+    // Deallocate all arrays
+    //
+    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t1,t2,diag,c,d,cf,d"
+                        "f)");
+    //
     // End of Cunhr_col01
     //
}
