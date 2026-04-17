--- a/mplapack/reference/Rtpttf.cpp
+++ b/mplapack/reference/Rtpttf.cpp
@@ -33,6 +33,10 @@
 //   Univ. of Colorado Denver
 //   NAG Ltd.
 
+#if defined(__INTEL_LLVM_COMPILER)
+#pragma clang optimize off
+#endif
+
 #include <mpblas.h>
 #include <mplapack.h>
 
