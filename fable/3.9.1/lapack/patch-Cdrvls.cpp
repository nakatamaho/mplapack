--- Cdrvls.cpp_	2026-01-20 06:55:53.286946187 +0900
+++ Cdrvls.cpp	2026-01-20 06:55:58.375065235 +0900
@@ -43,6 +43,8 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
+#include <memory>
+
 void Cdrvls(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, INTEGER const nnb, INTEGER *nbval, INTEGER *nxval, REAL const thresh, bool const tsterr, COMPLEX *a, COMPLEX *copya, COMPLEX *b, COMPLEX *copyb, COMPLEX *c, REAL *s, REAL *copys, INTEGER const nout) {
     common cmn;
     common_write write(cmn);
@@ -617,10 +619,6 @@
     //
     Alasvm(path, nout, nfail, nrun, nerrs);
     //
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(work)");
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(iwork)");
-    FEM_THROW_UNHANDLED("executable deallocate: deallocate(rwork)");
-    //
     // End of Cdrvls
     //
 }
