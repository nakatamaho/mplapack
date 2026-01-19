--- Cqrt04.cpp
+++ Cqrt04.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Cqrt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -93,7 +95,7 @@
     // Compute |I - Q'*Q| and store in RESULT(2)
     //
     Claset("Full", m, m, czero, one, r, m);
-    Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+    Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
     resid = Clansy("1", "Upper", m, r, m, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, m));
     //
@@ -184,9 +186,4 @@
     } else {
         result[6 - 1] = zero;
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
-    //
 }
