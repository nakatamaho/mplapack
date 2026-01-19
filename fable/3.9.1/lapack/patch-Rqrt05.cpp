--- Rqrt05.cpp
+++ Rqrt05.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Rqrt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -201,8 +203,4 @@
     } else {
         result[6 - 1] = zero;
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
 }
