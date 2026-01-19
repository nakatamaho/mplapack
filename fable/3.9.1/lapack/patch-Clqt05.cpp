--- Clqt05.cpp
+++ Clqt05.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Clqt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -97,8 +99,8 @@
     // Compute |L - A*Q*C| / |A| and store in RESULT(1)
     //
     Cgemm("N", "C", m, n2, n2, -one, a, m, q, n2, one, r, n2);
-    std::unique_ptr<COMPLEX[]> __rwork_storage(new COMPLEX[n2]);
-    COMPLEX *rwork = __rwork_storage.get();
+    std::unique_ptr<REAL[]> __rwork_storage(new REAL[n2]);
+    REAL *rwork = __rwork_storage.get();
     REAL anorm = Clange("1", m, n2, a, m, rwork);
     REAL resid = Clange("1", m, n2, r, n2, rwork);
     const REAL zero = 0.0;
@@ -111,7 +113,7 @@
     // Compute |I - Q*Q'| and store in RESULT(2)
     //
     Claset("Full", n2, n2, czero, one, r, n2);
-    Cherk("U", "N", n2, n2, dreal(-one), q, n2, dreal(one), r, n2);
+    Cherk("U", "N", n2, n2, (-one).real(), q, n2, one.real(), r, n2);
     resid = Clansy("1", "Upper", n2, r, n2, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, n2));
     //
@@ -204,8 +206,4 @@
     } else {
         result[6 - 1] = zero;
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
 }
