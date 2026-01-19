--- Rdrvls.cpp
+++ Rdrvls.cpp
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Rdrvls(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, INTEGER const nnb, INTEGER *nbval, INTEGER *nxval, REAL const thresh, bool const tsterr, REAL *a, REAL *copya, REAL *b, REAL *copyb, REAL *c, REAL *s, REAL *copys, INTEGER const nout) {
     common cmn;
     common_write write(cmn);
@@ -603,9 +605,6 @@
     //
     Alasvm(path, nout, nfail, nrun, nerrs);
     //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(work)");
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(iwork)");
-    //
     // End of Rdrvls
     //
 }
