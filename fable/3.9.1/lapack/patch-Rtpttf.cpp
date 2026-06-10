--- a/mplapack/reference/Rtpttf.cpp~	2026-02-05 08:12:25.350546391 +0900
+++ b/mplapack/reference/Rtpttf.cpp	2026-02-05 08:23:27.548010023 +0900
@@ -33,6 +33,10 @@
 //   Univ. of Colorado Denver
 //   NAG Ltd.
 
+#if defined(__INTEL_LLVM_COMPILER)
+#pragma clang optimize off
+#endif
+
 #include <mpblas.h>
 #include <mplapack.h>
 
