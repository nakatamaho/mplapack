--- Ctsqr01.cpp
+++ Ctsqr01.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Ctsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -56,14 +58,12 @@
     //
     REAL eps = Rlamch("Epsilon");
     INTEGER k = min(m, n);
-    INTEGER l = max(m, n, 1);
+    INTEGER l = max(m, n, (INTEGER)1);
     INTEGER mnb = max(mb, nb);
     INTEGER lwork = max((INTEGER)3, l) * mnb;
     //
     // Dynamically allocate local arrays
     //
-    // FABLE: ALLOCATE removed (RAII in C++)
-    //
     // Put random numbers into A and copy to AF
     //
     INTEGER j = 0;
@@ -101,8 +101,8 @@
     COMPLEX *q = __q_storage.get();
     std::unique_ptr<COMPLEX[]> __r_storage(new COMPLEX[m * l]);
     COMPLEX *r = __r_storage.get();
-    std::unique_ptr<COMPLEX[]> __rwork_storage(new COMPLEX[l]);
-    COMPLEX *rwork = __rwork_storage.get();
+    std::unique_ptr<REAL[]> __rwork_storage(new REAL[l]);
+    REAL *rwork = __rwork_storage.get();
     REAL anorm = 0.0;
     REAL resid = 0.0;
     const REAL zero = 0.0;
@@ -160,7 +160,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", m, m, czero, one, r, m);
-        Cherk("U", "C", m, m, dreal(-one), q, m, dreal(one), r, m);
+        Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
         resid = Clansy("1", "Upper", m, r, m, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, m));
         //
@@ -292,7 +292,7 @@
         // Compute |I - Q'*Q| and store in RESULT(2)
         //
         Claset("Full", n, n, czero, one, lq, l);
-        Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), lq, l);
+        Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), lq, l);
         resid = Clansy("1", "Upper", n, lq, l, rwork);
         result[2 - 1] = resid / (eps * max((INTEGER)1, n));
         //
@@ -377,9 +377,4 @@
         }
         //
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
-    //
 }
