--- Rlqt04.cpp
+++ Rlqt04.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Rlqt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
     common cmn;
     static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
@@ -183,9 +185,4 @@
     } else {
         result[6 - 1] = zero;
     }
-    //
-    // Deallocate all arrays
-    //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,l,rwork,work,t,c,d,cf,df)");
-    //
 }
