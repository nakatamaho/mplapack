--- a/TESTING/LIN/zlqt05.f~	2021-03-26 03:25:15.000000000 +0900
+++ b/TESTING/LIN/zlqt05.f	2026-01-28 07:33:44.525902792 +0900
@@ -93,8 +93,9 @@
 *     ..
 *     .. Local allocatable arrays
       COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:),
-     $  R(:,:), RWORK(:), WORK( : ), T(:,:),
+     $  R(:,:), WORK( : ), T(:,:),
      $  CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:)
+      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
 *
 *     .. Parameters ..
       DOUBLE PRECISION ZERO
