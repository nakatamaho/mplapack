https://github.com/Reference-LAPACK/lapack/pull/1203

diff --git a/SRC/chegvd.f b/SRC/chegvd.f
index 1a6ab4185..a20243ff9 100644
--- a/SRC/chegvd.f
+++ b/SRC/chegvd.f
@@ -319,7 +319,7 @@
 *
       IF( INFO.EQ.0 ) THEN
          WORK( 1 ) = SROUNDUP_LWORK(LOPT)
-         RWORK( 1 ) = REAL( LROPT )
+         RWORK( 1 ) = SROUNDUP_LWORK(LROPT)
          IWORK( 1 ) = LIOPT
 *
          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
@@ -357,9 +357,9 @@
       CALL CHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
      $             LRWORK,
      $             IWORK, LIWORK, INFO )
-      LOPT = INT( MAX( REAL( LOPT ), REAL( WORK( 1 ) ) ) )
-      LROPT = INT( MAX( REAL( LROPT ), REAL( RWORK( 1 ) ) ) )
-      LIOPT = INT( MAX( REAL( LIOPT ), REAL( IWORK( 1 ) ) ) )
+      LOPT = MAX( LOPT, INT( REAL( WORK( 1 ) ) ) )
+      LROPT = MAX( LROPT, INT( RWORK( 1 ) ) )
+      LIOPT = MAX( LIOPT, IWORK( 1 ) )
 *
       IF( WANTZ .AND. INFO.EQ.0 ) THEN
 *
@@ -396,7 +396,7 @@
       END IF
 *
       WORK( 1 ) = SROUNDUP_LWORK(LOPT)
-      RWORK( 1 ) = REAL( LROPT )
+      RWORK( 1 ) = SROUNDUP_LWORK(LROPT)
       IWORK( 1 ) = LIOPT
 *
       RETURN
diff --git a/SRC/chpgvd.f b/SRC/chpgvd.f
index ae55b4ee2..8a6634f2c 100644
--- a/SRC/chpgvd.f
+++ b/SRC/chpgvd.f
@@ -295,7 +295,7 @@
          END IF
 *
          WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
-         RWORK( 1 ) = REAL( LRWMIN )
+         RWORK( 1 ) = SROUNDUP_LWORK(LRWMIN)
          IWORK( 1 ) = LIWMIN
          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
             INFO = -11
@@ -331,9 +331,9 @@
       CALL CHPGST( ITYPE, UPLO, N, AP, BP, INFO )
       CALL CHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK,
      $             LRWORK, IWORK, LIWORK, INFO )
-      LWMIN = INT( MAX( REAL( LWMIN ), REAL( WORK( 1 ) ) ) )
-      LRWMIN = INT( MAX( REAL( LRWMIN ), REAL( RWORK( 1 ) ) ) )
-      LIWMIN = INT( MAX( REAL( LIWMIN ), REAL( IWORK( 1 ) ) ) )
+      LWMIN = MAX( LWMIN, INT( REAL( WORK( 1 ) ) ) )
+      LRWMIN = MAX( LRWMIN, INT( RWORK( 1 ) ) )
+      LIWMIN = MAX( LIWMIN, IWORK( 1 ) )
 *
       IF( WANTZ ) THEN
 *
@@ -377,7 +377,7 @@
       END IF
 *
       WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
-      RWORK( 1 ) = REAL( LRWMIN )
+      RWORK( 1 ) = SROUNDUP_LWORK(LRWMIN)
       IWORK( 1 ) = LIWMIN
       RETURN
 *
diff --git a/SRC/dspgvd.f b/SRC/dspgvd.f
index 7a0f6c073..b6bfb57d7 100644
--- a/SRC/dspgvd.f
+++ b/SRC/dspgvd.f
@@ -225,14 +225,15 @@
 *     ..
 *     .. External Functions ..
       LOGICAL            LSAME
-      EXTERNAL           LSAME
+      DOUBLE PRECISION   DROUNDUP_LWORK
+      EXTERNAL           LSAME, DROUNDUP_LWORK
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           DPPTRF, DSPEVD, DSPGST, DTPMV, DTPSV,
      $                   XERBLA
 *     ..
 *     .. Intrinsic Functions ..
-      INTRINSIC          DBLE, MAX
+      INTRINSIC          MAX
 *     ..
 *     .. Executable Statements ..
 *
@@ -268,7 +269,7 @@
                LWMIN = 2*N
             END IF
          END IF
-         WORK( 1 ) = LWMIN
+         WORK( 1 ) = DROUNDUP_LWORK(LWMIN)
          IWORK( 1 ) = LIWMIN
          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
             INFO = -11
@@ -302,8 +303,8 @@
       CALL DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
       CALL DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK,
      $             LIWORK, INFO )
-      LWMIN = INT( MAX( DBLE( LWMIN ), DBLE( WORK( 1 ) ) ) )
-      LIWMIN = INT( MAX( DBLE( LIWMIN ), DBLE( IWORK( 1 ) ) ) )
+      LWMIN = MAX( LWMIN, INT( WORK( 1 ) ) )
+      LIWMIN = MAX( LIWMIN, IWORK( 1 ) )
 *
       IF( WANTZ ) THEN
 *
@@ -346,7 +347,7 @@
          END IF
       END IF
 *
-      WORK( 1 ) = LWMIN
+      WORK( 1 ) = DROUNDUP_LWORK(LWMIN)
       IWORK( 1 ) = LIWMIN
 *
       RETURN
diff --git a/SRC/dsygvd.f b/SRC/dsygvd.f
index 6dbb777de..01181eb64 100644
--- a/SRC/dsygvd.f
+++ b/SRC/dsygvd.f
@@ -245,14 +245,15 @@
 *     ..
 *     .. External Functions ..
       LOGICAL            LSAME
-      EXTERNAL           LSAME
+      DOUBLE PRECISION   DROUNDUP_LWORK
+      EXTERNAL           LSAME, DROUNDUP_LWORK
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           DPOTRF, DSYEVD, DSYGST, DTRMM, DTRSM,
      $                   XERBLA
 *     ..
 *     .. Intrinsic Functions ..
-      INTRINSIC          DBLE, MAX
+      INTRINSIC          MAX
 *     ..
 *     .. Executable Statements ..
 *
@@ -290,7 +291,7 @@
       END IF
 *
       IF( INFO.EQ.0 ) THEN
-         WORK( 1 ) = LOPT
+         WORK( 1 ) = DROUNDUP_LWORK(LOPT)
          IWORK( 1 ) = LIOPT
 *
          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
@@ -326,8 +327,8 @@
       CALL DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
      $             LIWORK,
      $             INFO )
-      LOPT = INT( MAX( DBLE( LOPT ), DBLE( WORK( 1 ) ) ) )
-      LIOPT = INT( MAX( DBLE( LIOPT ), DBLE( IWORK( 1 ) ) ) )
+      LOPT = MAX( LOPT, INT( WORK( 1 ) ) )
+      LIOPT = MAX( LIOPT, IWORK( 1 ) )
 *
       IF( WANTZ .AND. INFO.EQ.0 ) THEN
 *
@@ -363,7 +364,7 @@
          END IF
       END IF
 *
-      WORK( 1 ) = LOPT
+      WORK( 1 ) = DROUNDUP_LWORK(LOPT)
       IWORK( 1 ) = LIOPT
 *
       RETURN
diff --git a/SRC/sspgvd.f b/SRC/sspgvd.f
index 91067e418..b4780aa73 100644
--- a/SRC/sspgvd.f
+++ b/SRC/sspgvd.f
@@ -233,7 +233,7 @@
      $                   XERBLA
 *     ..
 *     .. Intrinsic Functions ..
-      INTRINSIC          MAX, REAL
+      INTRINSIC          MAX
 *     ..
 *     .. Executable Statements ..
 *
@@ -303,8 +303,8 @@
       CALL SSPGST( ITYPE, UPLO, N, AP, BP, INFO )
       CALL SSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK,
      $             LIWORK, INFO )
-      LWMIN = INT( MAX( REAL( LWMIN ), REAL( WORK( 1 ) ) ) )
-      LIWMIN = INT( MAX( REAL( LIWMIN ), REAL( IWORK( 1 ) ) ) )
+      LWMIN = MAX( LWMIN, INT( WORK( 1 ) ) )
+      LIWMIN = MAX( LIWMIN, IWORK( 1 ) )
 *
       IF( WANTZ ) THEN
 *
diff --git a/SRC/ssygvd.f b/SRC/ssygvd.f
index aacf9d3ad..2ca694258 100644
--- a/SRC/ssygvd.f
+++ b/SRC/ssygvd.f
@@ -253,7 +253,7 @@
      $                   XERBLA
 *     ..
 *     .. Intrinsic Functions ..
-      INTRINSIC          MAX, REAL
+      INTRINSIC          MAX
 *     ..
 *     .. Executable Statements ..
 *
@@ -327,8 +327,8 @@
       CALL SSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
      $             LIWORK,
      $             INFO )
-      LOPT = INT( MAX( REAL( LOPT ), REAL( WORK( 1 ) ) ) )
-      LIOPT = INT( MAX( REAL( LIOPT ), REAL( IWORK( 1 ) ) ) )
+      LOPT = MAX( LOPT, INT( WORK( 1 ) ) )
+      LIOPT = MAX( LIOPT, IWORK( 1 ) )
 *
       IF( WANTZ .AND. INFO.EQ.0 ) THEN
 *
diff --git a/SRC/zhegvd.f b/SRC/zhegvd.f
index 208735b0e..bdcd681a9 100644
--- a/SRC/zhegvd.f
+++ b/SRC/zhegvd.f
@@ -268,7 +268,8 @@
 *     ..
 *     .. External Functions ..
       LOGICAL            LSAME
-      EXTERNAL           LSAME
+      DOUBLE PRECISION   DROUNDUP_LWORK
+      EXTERNAL           LSAME, DROUNDUP_LWORK
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           XERBLA, ZHEEVD, ZHEGST, ZPOTRF, ZTRMM,
@@ -317,8 +318,8 @@
       END IF
 *
       IF( INFO.EQ.0 ) THEN
-         WORK( 1 ) = LOPT
-         RWORK( 1 ) = REAL( LROPT )
+         WORK( 1 ) = DROUNDUP_LWORK(LOPT)
+         RWORK( 1 ) = DROUNDUP_LWORK(LROPT)
          IWORK( 1 ) = LIOPT
 *
          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
@@ -356,9 +357,9 @@
       CALL ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
      $             LRWORK,
      $             IWORK, LIWORK, INFO )
-      LOPT = INT( MAX( DBLE( LOPT ), DBLE( WORK( 1 ) ) ) )
-      LROPT = INT( MAX( DBLE( LROPT ), DBLE( RWORK( 1 ) ) ) )
-      LIOPT = INT( MAX( DBLE( LIOPT ), DBLE( IWORK( 1 ) ) ) )
+      LOPT = MAX( LOPT, INT( DBLE( WORK( 1 ) ) ) )
+      LROPT = MAX( LROPT, INT( RWORK( 1 ) ) )
+      LIOPT = MAX( LIOPT, IWORK( 1 ) )
 *
       IF( WANTZ .AND. INFO.EQ.0 ) THEN
 *
@@ -394,8 +395,8 @@
          END IF
       END IF
 *
-      WORK( 1 ) = LOPT
-      RWORK( 1 ) = REAL( LROPT )
+      WORK( 1 ) = DROUNDUP_LWORK(LOPT)
+      RWORK( 1 ) = DROUNDUP_LWORK(LROPT)
       IWORK( 1 ) = LIOPT
 *
       RETURN
diff --git a/SRC/zhpgvd.f b/SRC/zhpgvd.f
index 1f1041ffc..78414fa5f 100644
--- a/SRC/zhpgvd.f
+++ b/SRC/zhpgvd.f
@@ -246,7 +246,8 @@
 *     ..
 *     .. External Functions ..
       LOGICAL            LSAME
-      EXTERNAL           LSAME
+      DOUBLE PRECISION   DROUNDUP_LWORK
+      EXTERNAL           LSAME, DROUNDUP_LWORK
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           XERBLA, ZHPEVD, ZHPGST, ZPPTRF, ZTPMV,
@@ -293,8 +294,8 @@
             END IF
          END IF
 *
-         WORK( 1 ) = LWMIN
-         RWORK( 1 ) = REAL( LRWMIN )
+         WORK( 1 ) = DROUNDUP_LWORK(LWMIN)
+         RWORK( 1 ) = DROUNDUP_LWORK(LRWMIN)
          IWORK( 1 ) = LIWMIN
          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
             INFO = -11
@@ -330,9 +331,9 @@
       CALL ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
       CALL ZHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK,
      $             LRWORK, IWORK, LIWORK, INFO )
-      LWMIN = INT( MAX( DBLE( LWMIN ), DBLE( WORK( 1 ) ) ) )
-      LRWMIN = INT( MAX( DBLE( LRWMIN ), DBLE( RWORK( 1 ) ) ) )
-      LIWMIN = INT( MAX( DBLE( LIWMIN ), DBLE( IWORK( 1 ) ) ) )
+      LWMIN = MAX( LWMIN, INT( DBLE( WORK( 1 ) ) ) )
+      LRWMIN = MAX( LRWMIN, INT( RWORK( 1 ) ) )
+      LIWMIN = MAX( LIWMIN, IWORK( 1 ) )
 *
       IF( WANTZ ) THEN
 *
@@ -375,8 +376,8 @@
          END IF
       END IF
 *
-      WORK( 1 ) = LWMIN
-      RWORK( 1 ) = REAL( LRWMIN )
+      WORK( 1 ) = DROUNDUP_LWORK(LWMIN)
+      RWORK( 1 ) = DROUNDUP_LWORK(LRWMIN)
       IWORK( 1 ) = LIWMIN
       RETURN
 *
