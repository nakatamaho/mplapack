--- Clqt04.cpp
+++ Clqt04.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Clqt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -79,8 +81,8 @@
     // Compute |L - A*Q'| / |A| and store in RESULT(1)
     //
     Cgemm("N", "C", m, n, n, -one, a, m, q, n, one, l, ll);
-    std::unique_ptr<COMPLEX[]> __rwork_storage(new COMPLEX[ll]);
-    COMPLEX *rwork = __rwork_storage.get();
+    std::unique_ptr<REAL[]> __rwork_storage(new REAL[ll]);
+    REAL *rwork = __rwork_storage.get();
     REAL anorm = Clange("1", m, n, a, m, rwork);
     REAL resid = Clange("1", m, n, l, ll, rwork);
     const REAL zero = 0.0;
@@ -93,7 +95,7 @@
     // Compute |I - Q'*Q| and store in RESULT(2)
     //
     Claset("Full", n, n, czero, one, l, ll);
-    Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), l, ll);
+    Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), l, ll);
     resid = Clansy("1", "Upper", n, l, ll, rwork);
     result[2 - 1] = resid / (eps * max((INTEGER)1, n));
     //
@@ -184,9 +186,4 @@
     } else {
         result[6 - 1] = zero;
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,l,rwork,work,t,c,d,cf,df)");
-    //
 }
