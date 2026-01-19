--- Cqrt05.cpp
+++ Cqrt05.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Cqrt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -111,7 +113,7 @@
     // Compute |I - Q'*Q| and store in RESULT(2)
     //
     Claset("Full", m2, m2, czero, one, r, m2);
-    Cherk("U", "C", m2, m2, dreal(-one), q, m2, dreal(one), r, m2);
+    Cherk("U", "C", m2, m2, -one.real(), q, m2, one.real(), r, m2);
     resid = Clansy("1", "Upper", m2, r, m2, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m2));
     //
@@ -202,8 +204,4 @@
     } else {
         result[6 - 1] = zero;
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
 }
