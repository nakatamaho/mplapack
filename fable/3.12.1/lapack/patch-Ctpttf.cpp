--- a/mplapack/reference/Ctpttf.cpp~	2026-02-05 08:13:17.412401988 +0900
+++ b/mplapack/reference/Ctpttf.cpp	2026-02-05 08:18:09.038366746 +0900
@@ -33,6 +33,10 @@
 //   Univ. of Colorado Denver
 //   NAG Ltd.
 
+#if defined(__INTEL_LLVM_COMPILER)
+#pragma clang optimize off
+#endif
+
 #include <mpblas.h>
 #include <mplapack.h>
 
