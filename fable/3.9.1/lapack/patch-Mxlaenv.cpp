--- Mxlaenv.cpp
+++ Mxlaenv.cpp
@@ -44,8 +44,6 @@
 #include <mplapack_lin.h>
 
 void Mxlaenv(INTEGER const ispec, INTEGER const nvalue) {
-    common cmn;
-    arr_ref<int> iparms(cmn.iparms, dimension(100));
     //
     if (ispec >= 1 && ispec <= 9) {
         iparms[ispec - 1] = nvalue;
