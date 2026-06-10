--- a/TESTING/LIN/ztsqr01.f~	2021-03-26 03:25:15.000000000 +0900
+++ b/TESTING/LIN/ztsqr01.f	2026-01-28 07:37:57.411322337 +0900
@@ -96,8 +96,9 @@
 *     ..
 *     .. Local allocatable arrays
       COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:),
-     $  R(:,:), RWORK(:), WORK( : ), T(:),
+     $  R(:,:), WORK( : ), T(:),
      $  CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:), LQ(:,:)
+      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
 *
 *     .. Parameters ..
       DOUBLE PRECISION ZERO
