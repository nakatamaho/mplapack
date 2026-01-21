--- Clahilb.cpp_	2026-01-21 19:59:05.683503438 +0900
+++ Clahilb.cpp	2026-01-21 19:59:09.752575846 +0900
@@ -43,7 +43,6 @@
 #include <mplapack_matgen.h>
 
 void Clahilb(INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, REAL *work, INTEGER &info, fem::str_cref path) {
-    common cmn;
     static COMPLEX d1[8] = {COMPLEX(-1.0, 0.0), COMPLEX(0.0, 1.0), COMPLEX(-1.0, -1.0), COMPLEX(0.0, -1.0), COMPLEX(1.0, 0.0), COMPLEX(-1.0, 1.0), COMPLEX(1.0, 1.0), COMPLEX(1.0, -1.0)};
     static COMPLEX d2[8] = {COMPLEX(-1.0, 0.0), COMPLEX(0.0, -1.0), COMPLEX(-1.0, 1.0), COMPLEX(0.0, 1.0), COMPLEX(1.0, 0.0), COMPLEX(-1.0, -1.0), COMPLEX(1.0, -1.0), COMPLEX(1.0, 1.0)};
     static COMPLEX invd1[8] = {COMPLEX(-1.0, 0.0), COMPLEX(0.0, -1.0), COMPLEX(-0.5, 0.5), COMPLEX(0.0, 1.0), COMPLEX(1.0, 0.0), COMPLEX(-0.5, -0.5), COMPLEX(0.5, -0.5), COMPLEX(0.5, 0.5)};
