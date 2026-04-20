--- a/TESTING/LIN/zlqt04.f_	2026-01-28 07:11:40.162066743 +0900
+++ b/TESTING/LIN/zlqt04.f	2026-01-28 07:12:42.875336664 +0900
@@ -86,8 +86,9 @@
 *     ..
 *     .. Local allocatable arrays
       COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:),
-     $  L(:,:), RWORK(:), WORK( : ), T(:,:),
+     $  L(:,:), WORK( : ), T(:,:),
      $  CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:)
+      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
 *
 *     .. Parameters ..
       DOUBLE PRECISION ZERO
