--- Rchkaa.cpp_	2026-01-20 13:07:14.185065602 +0900
+++ Rchkaa.cpp	2026-01-20 13:08:31.590863365 +0900
@@ -43,8 +43,7 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_dchkaa(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+void Rchkaa(void) {
     common_read read(cmn);
     common_write write(cmn);
     static fem::str<10> intstr = "0123456789";
@@ -934,4 +933,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkaa); }
+int main(int argc, char const *argv[]) { return Rchkaa(); }
