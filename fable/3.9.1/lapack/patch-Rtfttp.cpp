--- a/mplapack/reference/Rtfttp.cpp~	2026-02-05 08:12:25.027541118 +0900
+++ b/mplapack/reference/Rtfttp.cpp	2026-02-05 08:22:01.347466902 +0900
@@ -33,6 +33,10 @@
 //   Univ. of Colorado Denver
 //   NAG Ltd.
 
+#if defined(__INTEL_LLVM_COMPILER)
+#pragma clang optimize off
+#endif
+
 #include <mpblas.h>
 #include <mplapack.h>
 
