dnl
dnl Configure for MPACK
dnl

AC_INIT(mpack, 0.8.0)
AC_CONFIG_HEADER(include/config.h) 

AC_CONFIG_FILES([ Makefile mpack.pc ])

AC_PROG_CC
AC_PROG_CXX

AC_PROG_LIBTOOL
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_MACRO_DIR([m4])

dnl OpenMP Compiler flags.
AC_OPENMP
AC_LANG_PUSH([C++])
AC_OPENMP
AC_LANG_POP

dnl Checking host name
AC_CANONICAL_HOST
IS_MINGW=0
case "$host" in
    *-mingw32*)
        IS_MINGW=1
        ;;
    *)
esac
AM_CONDITIONAL(IS_MINGW, test x$IS_MINGW = x1)
IS_MACOSX=0
case "$host" in
    *apple-darwin*)
        IS_MACOSX=1
        ;;
    *)
esac
AM_CONDITIONAL(IS_MACOSX, test x$IS_MACOSX = x1)

IS_LINUX=0
case "$host" in
    *linux*)
        IS_LINUX=1
        ;;
    *)
esac
DYLD=
AM_CONDITIONAL(IS_LINUX, test x$IS_LINUX = x1)
if test x"$IS_LINUX" = x"1"; then
DYLD=-ldl
fi
AC_SUBST(DYLD)

dnl dummy for Intel Xeon Phi (aka MIC)
IS_MIC=0
AM_CONDITIONAL(IS_MIC, test x$IS_MIC = x1)

dnl enable debug; compare with reference LAPACK and BLAS. 
dnl when enable this flag, we perform full build;
dnl all mp versions of mlapack and mblas, etc.
AC_ARG_ENABLE(debug,
[  --enable-debug    Enable debug],
[case "${enableval}" in
    yes) enable_debug=yes ;;
    no)  enable_debug=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
        esac],[enable_debug=no])
AM_CONDITIONAL(ENABLE_DEBUG, test x$enable_debug = xyes)

AC_ARG_ENABLE(optimization,
[  --enable-optimization    Enable compiler level optimization],
[case "${enableval}" in
    yes) enable_optimization=yes ;;
    no)  enable_optimization=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-optimization) ;;
        esac],[enable_optimization=yes])

AC_ARG_ENABLE(cuda,
[  --enable-cuda    Enable compiler level cuda],
[case "${enableval}" in
    yes) enable_cuda=yes ;;
    no)  enable_cuda=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-cuda) ;;
        esac],[enable_cuda=no])
AM_CONDITIONAL(ENABLE_CUDA, test x$enable_cuda = xyes)

AC_ARG_WITH(cudatoolkithome,
[  --with-cudatoolkithome     Specify path to CUDA Toolkit (default: /usr/local/cuda)],
[ CUDATOOLKITHOME="$withval"], 
[ CUDATOOLKITHOME=/usr/local/cuda]
)

ARCH=`uname -m`
if test $ARCH = "x86_64"; then
  SUFFIX="64"
else
  SUFFIX=""
fi

CUDA_CFLAGS="-I$CUDATOOLKITHOME/include"
CUDA_LIBS="-L$CUDATOOLKITHOME/lib$SUFFIX"
CUDA_LDFLAGS="-L$CUDATOOLKITHOME/lib$SUFFIX"
NVCC="$CUDATOOLKITHOME/bin/nvcc"

IS_CUDA_VERSION_MAJOR_3="FALSE"
dnl we should write check part for CUDA toolkit
if test x$enable_cuda = x"yes"; then
CUDA_VERSION=`"$NVCC" --version | grep "Cuda compilation tools"  | awk '{print $5}' | sed 's/,//g'`
CUDA_VERSION_MAJOR=`"$NVCC" --version | grep "Cuda compilation tools"  | awk '{print $5}' | sed 's/,//g' | sed 's/\..$//g'`
AC_MSG_RESULT([CUDA version is $CUDA_VERSION])
AC_MSG_RESULT([CUDA major version is $CUDA_VERSION_MAJOR])
fi

AM_CONDITIONAL(IS_CUDA_VERSION_MAJOR_3, test x"$CUDA_VERSION_MAJOR" = x"3")

AC_SUBST(CUDATOOLKITHOME)
AC_SUBST(CUDA_LDFLAGS)
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(CUDA_VERSION)
AC_SUBST(CUDA_VERSION_MAJOR)
AC_SUBST(IS_CUDA_VERSION_MAJOR_3)
AC_SUBST(NVCC)

dnl enable FORTRAN warpper
AC_ARG_ENABLE(fortranwrapper,
[  --enable-fortranwrapper    Enable fortranwrapper],
[case "${enableval}" in
    yes) enable_fortranwrapper=yes ;;
    no)  enable_fortranwrapper=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-fortranwrapper) ;;
        esac],[enable_fortranwrapper=no])
AM_CONDITIONAL(ENABLE_FORTRANWRAPPER, test x$enable_fortranwrapper = xyes)

AC_ARG_ENABLE(wrapper,
[  --enable-wrapper    Enable wrapper for BLAS and LAPACK. xBLAS like approach],
[case "${enableval}" in
    yes) with_wrapper=yes ;;
    no)  with_wrapper=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-wrapper) ;;
        esac],[with_wrapper=no])
AM_CONDITIONAL(ENABLE_WRAPPER, test x$with_wrapper = xyes)

dnl enable GMP version of MPACK.
AC_ARG_ENABLE(gmp,
[  --enable-gmp    Enable GMP version],
[case "${enableval}" in
    yes) enable_gmp=yes ;;
    no)  enable_gmp=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-gmp) ;;
        esac],[enable_gmp=yes])
AM_CONDITIONAL(ENABLE_GMP, test x$enable_gmp = xyes)

dnl enable MPFR version of MPACK.
AC_ARG_ENABLE(mpfr,
[  --enable-mpfr    Enable MPFR version],
[case "${enableval}" in
    yes) enable_mpfr=yes ;;
    no)  enable_mpfr=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpfr) ;;
        esac],[enable_mpfr=no])

dnl We use MPFR version as reference
dnl so if debug is turned on, 
if test x$enable_debug = x"yes"; then
    enable_mpfr=yes
fi
AM_CONDITIONAL(ENABLE_MPFR, test x$enable_mpfr = xyes)

dnl enable QD version of MPACK.
AC_ARG_ENABLE(qd,
[  --enable-qd    Enable QD version],
[case "${enableval}" in
    yes) enable_qd=yes ;;
    no)  enable_qd=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-qd) ;;
        esac],[enable_qd=no])
AM_CONDITIONAL(ENABLE_QD, test x$enable_qd = xyes)

AC_ARG_WITH(qdlib-flags,
[  --with-qdlib-flags      passing flag to internal qd library
],,with_qdlib_flags="")
QDLIB_FLAGS=$with_qdlib_flags;
AC_SUBST(QDLIB_FLAGS)

dnl enable DD version of MPACK.
AC_ARG_ENABLE(dd,
[  --enable-dd    Enable DD version],
[case "${enableval}" in
    yes) enable_dd=yes ;;
    no)  enable_dd=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-dd) ;;
        esac],[ 
if test x"$IS_MINGW" = x"1"; then
enable_dd=no
else
enable_dd=yes
fi
])
AM_CONDITIONAL(ENABLE_DD, test x$enable_dd = xyes)

dnl enable double
AC_ARG_ENABLE(double,
[  --enable-double    Enable double version],
[case "${enableval}" in
    yes) enable_double=yes ;;
    no)  enable_double=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-double) ;;
        esac],[enable_double=no])
AM_CONDITIONAL(ENABLE_DOUBLE, test x$enable_double = xyes)

AC_ARG_ENABLE(__float128,
[  --enable-__float128    Enable __float128 version using non standard extension of gcc and icc],
[case "${enableval}" in
    yes) enable___float128=yes ;;
    no)  enable___float128=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-__float128) ;;
        esac],[enable___float128=no])

AC_ARG_WITH(gccname,
[  --with-gccname     Specify path to gcc for libquadmath (only for Intel Composer) ],
[ GCCNAME="$withval"]
)

if test x$enable___float128 = x"yes"; then

__ICC=`basename $CC`
if test x$__ICC = x"icc"; then
AC_MSG_RESULT([Note: your compiler is Intel Composer, so you need -lquadmath, thus -gcc-name])
if test x$GCCNAME != x"" ; then
    _GCC_FLOAT128_FLAGS="-gcc-name=$GCCNAME"
else
__GCC=`which gcc`
    _GCC_FLOAT128_FLAGS="-gcc-name=$__GCC"
fi
fi

SAVE_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $_GCC_FLOAT128_FLAGS"
AC_MSG_CHECKING([whether __float128 is supported])
AC_TRY_COMPILE([
#include <quadmath.h>
],[__float128 a = 0.0Q; return 0;], have___float128=yes, have___float128=no)

if test "$have___float128" = "no"; then
    AC_MSG_ERROR([__float128 is not supported on this platform])
else
    AC_MSG_RESULT([__float128 is supported])
fi
fi
CFLAGS="$SAVE_CFLAGS"

AM_CONDITIONAL(ENABLE___FLOAT128, test x$enable___float128 = xyes)

dnl GMP library start
AC_ARG_WITH(system-gmp,
[  --with-system-gmp      Use GMP library already on system
],,with_system_gmp="no")
GMP_INCLUDEDIR=.
GMP_LIBDIR=.
if test x$enable_gmp = x"yes" -o x$enable_mpfr = x"yes" -o x$enable_debug = x"yes"; then
    if test x$with_system_gmp = x"yes"; then
        HAVE_GMP=""
        AC_CHECK_HEADER(gmp.h, HAVE_GMP_HEADER=yes)
        if test x$HAVE_GMP_HEADER = xyes; then
            save_LDFLAGS="$LDFLAGS"
            LDFLAGS="$LDFLAGS -L$GMP_LIBDIR"
            AC_CHECK_LIB(gmp, __gmpz_init, [GMP_LIBS="-lgmpxx -lgmp"; HAVE_GMP=yes])
            LDFLAGS=$save_LDFLAGS
        fi
        if test x$HAVE_GMP != xyes; then
            AC_MSG_RESULT([No GMP library with C++ wrapper found; using internal one instead])
        else
            AC_MSG_RESULT([GMP library found: using GMP library comes with the system])
            SYSTEM_GMP=yes
        fi
    else
        AC_MSG_RESULT([Using GMP library comes with MPACK])
        GMP_INCLUDEDIR=external/i/GMP/include/
        GMP_LIBDIR=external/i/GMP/lib/
    fi
fi
AM_CONDITIONAL(SYSTEM_GMP, test x$SYSTEM_GMP = xyes)
AC_SUBST(GMP_INCLUDEDIR)
AC_SUBST(GMP_LIBDIR)
dnl GMP library end

dnl MPFR library start
AC_ARG_WITH(system-mpfr,
[  --with-system-mpfr      Use MPFR library already on system
],,with_system_mpfr="no")
MPFR_INCLUDEDIR=.
MPFR_LIBDIR=.
if test x$enable_mpfr = x"yes" -o x$enable_debug = x"yes"; then
    if test x$with_system_mpfr = x"yes"; then
        HAVE_MPFR=""
        AC_CHECK_HEADER(mpfr.h, HAVE_MPFR_HEADER=yes)
        if test x$HAVE_MPFR_HEADER = xyes; then
            save_LDFLAGS="$LDFLAGS"
            LDFLAGS="$LDFLAGS -L$MPFR_LIBDIR"
            AC_CHECK_LIB(mpfr, mpfr_init, [MPFR_LIBS="-lmpfr -lgmp"; HAVE_MPFR=yes])
            LDFLAGS=$save_LDFLAGS
        fi
        if test x$HAVE_MPFR != xyes; then
            AC_MSG_RESULT([No MPFR library found; using internal one instead])
        else
            AC_MSG_RESULT([MPFR found: using MPFR comes with the system])
            SYSTEM_MPFR=yes
        fi
    else
        AC_MSG_RESULT([Using MPFR comes with MPACK])
        MPFR_INCLUDEDIR=external/i/MPFR/include/
        MPFR_LIBDIR=external/i/MPFR/lib/
    fi
fi
AM_CONDITIONAL(SYSTEM_MPFR, test x$SYSTEM_MPFR = xyes)
AC_SUBST(MPFR_INCLUDEDIR)
AC_SUBST(MPFR_LIBDIR)
dnl MPFR library end

dnl MPC library start
AC_ARG_WITH(system-mpc,
[  --with-system-mpc      Use MPC library already on system
],,with_system_mpc="no")
MPC_INCLUDEDIR=.
MPC_LIBDIR=.
if test x$enable_mpfr = x"yes" -o x$enable_debug = x"yes"; then
    if test x$with_system_mpc = x"yes"; then
        HAVE_MPC=""
        AC_CHECK_HEADER(mpc.h, HAVE_MPC_HEADER=yes)
        if test x$HAVE_MPC_HEADER = xyes; then
            save_LDFLAGS="$LDFLAGS"
            LDFLAGS="$LDFLAGS -L$MPC_LIBDIR"
            AC_CHECK_LIB(mpc, mpc_init2, [MPC_LIBS="-lmpc -lmpfr -lgmp"; HAVE_MPC=yes])
            LDFLAGS=$save_LDFLAGS
        fi
        if test x$HAVE_MPC != xyes; then
            AC_MSG_RESULT([No MPC library found; using internal one instead])
        else
            AC_MSG_RESULT([MPC library found: using MPC library comes with the system])
            SYSTEM_MPC=yes
	    MPC_INCLUDEDIR=.
	    MPC_LIBDIR=.
        fi
    else
	AC_MSG_RESULT([Using MPC library comes with MPACK])
        MPC_INCLUDEDIR=external/i/MPC/include/
        MPC_LIBDIR=external/i/MPC/lib/
    fi
fi
AM_CONDITIONAL(SYSTEM_MPC, test x$SYSTEM_MPC = xyes)
AC_SUBST(MPC_INCLUDEDIR)
AC_SUBST(MPC_LIBDIR)
dnl MPC library end

dnl QD library start
AC_ARG_WITH(system-qd,
[  --with-system-qd      Use QD library already on system
],,with_system_qd="no")
QD_INCLUDEDIR=.
QD_LIBDIR=.
if test x$enable_qd = x"yes" -o x$enable_dd = x"yes" -o x$enable_debug = x"yes"; then
    if test x$with_system_qd = x"yes"; then
        HAVE_QD=""
        AC_LANG_PUSH([C++])
        AC_CHECK_HEADER(qd/qd_real.h, HAVE_QD_HEADER=yes)
        if test x$HAVE_QD_HEADER = xyes; then
            save_LDFLAGS="$LDFLAGS"
            LDFLAGS="$LDFLAGS -L$QD_LIBDIR"
            AC_CHECK_LIB(qd, c_qd_sqrt, [QD_LIBS="-lqd"; HAVE_QD=yes])
            LDFLAGS=$save_LDFLAGS
        fi
        if test x$HAVE_QD != xyes; then
              AC_MSG_RESULT([No QD library found; using internal one instead])
        else
              AC_MSG_RESULT([QD library found: using QD library comes with the system])
              SYSTEM_QD=yes
	      QD_INCLUDEDIR=.
	      QD_LIBDIR=.
        fi
        AC_LANG_POP([C++])
    else
        AC_MSG_RESULT([Using QD library comes with MPACK])
        QD_INCLUDEDIR=external/i/QD/include/
        QD_LIBDIR=external/i/QD/lib/
    fi   
fi
AM_CONDITIONAL(SYSTEM_QD, test x$SYSTEM_QD = xyes)
AC_SUBST(QD_INCLUDEDIR)
AC_SUBST(QD_LIBDIR)
dnl QD library end

need_internal_gmplib=no
if test x$enable_gmp = x"yes" && test x$SYSTEM_GMP != x"yes"; then
    need_internal_gmplib=yes
fi
AM_CONDITIONAL(NEED_INTERNAL_GMPLIB, test x$need_internal_gmplib = xyes)

need_internal_qdlib=no
if test x$enable_qd = x"yes" || test x$enable_dd = x"yes"; then
    if test x$SYSTEM_QD != x"yes"; then
       need_internal_qdlib=yes
    fi
fi
AM_CONDITIONAL(NEED_INTERNAL_QDLIB, test x$need_internal_qdlib = xyes)

need_internal_mpfrlib=no
if test x$enable_mpfr = x"yes" && test x$SYSTEM_MPFR != x"yes"; then
    need_internal_mpfrlib=yes
fi
AM_CONDITIONAL(NEED_INTERNAL_MPFRLIB, test x$need_internal_mpfrlib = xyes)

need_internal_mpclib=no
if test x$enable_mpfr = x"yes" && test x$SYSTEM_MPC != x"yes"; then
    need_internal_mpclib=yes
fi
AM_CONDITIONAL(NEED_INTERNAL_MPCLIB, test x$need_internal_mpclib = xyes)

dnl print status...start
if test x$enable_gmp = x"yes"; then
AC_MSG_RESULT([Enable GMP version])
else
AC_MSG_RESULT([Disable GMP version])
fi
if test x$enable_dd = x"yes"; then
AC_MSG_RESULT([Enable DD version])
else
AC_MSG_RESULT([Disable DD version])
fi

if test x$enable_cuda = x"yes"; then
AC_MSG_RESULT([Enable CUDA version : the home of CUDA Toolkit : $CUDATOOLKITHOME])
else
AC_MSG_RESULT([Disable CUDA version])
fi

if test x"$IS_MINGW" = x"1" && test x$enable_dd = x"yes"; then
#it builds fine but we have to explicitly call mpack_initialize_qd(void);
#when we start dd calculation and call mpack_finalize_qd(void)
#at the end of calculation
    AC_MSG_ERROR([dd is not supported on this platform; constructor issue])
fi

if test x$enable_qd = x"yes"; then
AC_MSG_RESULT([Enable QD version])
else
AC_MSG_RESULT([Disable QD version])
fi

if test x"$IS_MINGW" = x"1" && test x$enable_qd = x"yes"; then
#it builds fine but we have to explicitly call mpack_initialize_qd(void);
#when we start dd calculation and call mpack_finalize_qd(void)
#at the end of calculation
    AC_MSG_ERROR([qd is not supported on this platform; constructor issue])
fi

if test x$enable_mpfr = x"yes"; then
AC_MSG_RESULT([Enable MPFR version])
else
AC_MSG_RESULT([Disable MPFR version])
fi
if test x$enable_double = x"yes"; then
AC_MSG_RESULT([Enable double version])
else
AC_MSG_RESULT([Disable double version])
fi
if test x$enable___float128 = x"yes"; then
AC_MSG_RESULT([Enable __float128 version])
else
AC_MSG_RESULT([Disable __float128 version])
fi
if test x$enable_debug = x"yes"; then
AC_MSG_RESULT([Enable debug])
else
AC_MSG_RESULT([Disable debug])
fi
if test x$enable_fortranwrapper = x"yes"; then
AC_MSG_RESULT([Enable fortran wrapper])
else
AC_MSG_RESULT([Disable fortran wrapper])
fi

dnl print status...end

AC_ARG_WITH(external-blas,
[  --with-external-blas           Specify location of external BLAS e.g. "-L/usr/local/lib -lblas". It is used for benchmark only.
],,)

AC_ARG_WITH(external-lapack,
[  --with-external-lapack         Specify how to link against external LAPACK e.g. "-L/usr/local/lib -llapack". It is used for benchmark only.
],,)

dnl First we check Fotran 90 compiler stuffs.

if test x"$enable_debug" = x"yes" || test x"$enable_fortranwrapper" = x"yes"; then
    AC_PROG_FC
    AC_FC_LIBRARY_LDFLAGS
    AC_FC_WRAPPERS
    AC_SUBST(FCLIBS)
    AC_LANG_PUSH([Fortran])
    AC_OPENMP
    AC_LANG_POP([Fortran])
fi

dnl For benchmarking, we can set explicitly BLAS and LAPACK
dnl libraries, like Intel MKL, GotoBLAS2, etc.
dnl very specific case: veclib is shipped for MacOSX. it should be used without fortran compiler
dnl this library should be the default for benchmarking (actually it is ATLAS).
dnl TODO: MacOSX support.
if test x"$with_external_blas" != x""; then
    HAVE_EXTERNAL_BLAS=""
    BLAS=$with_external_blas;
    save_LIBS="$LIBS"
    AC_FC_FUNC(dgemm)
    LIBS="$LIBS $BLAS $FCLIBS"
    AC_CHECK_FUNC($dgemm, [HAVE_EXTERNAL_BLAS=yes])
    LIBS="$save_LIBS"
    AC_SUBST(BLAS)
    if test x"$HAVE_EXTERNAL_BLAS" != xyes; then
        AC_MSG_ERROR([No external BLAS found])
    fi
fi
if test x"$with_external_lapack" != x""; then
    dnl Check whether we can call LAPACK subroutines
    HAVE_EXTERNAL_LAPACK=""
    LAPACK=$with_external_lapack;
    dnl default blas name is libblas, thus we pass -lblas.
    if test x"$with_external_blas" = x; then
        BLAS="-lblas"
    fi
    save_LIBS="$LIBS"
    AC_FC_FUNC(dsyev)
    LIBS="$LIBS $LAPACK $BLAS $FCLIBS"
    AC_CHECK_FUNC($dsyev, [HAVE_EXTERNAL_LAPACK=yes])
    LIBS="$save_LIBS"
    AC_SUBST(LAPACK)
    if test x"$HAVE_EXTERNAL_LAPACK" != xyes; then
       AC_MSG_ERROR([No external LAPACK found])
    fi
fi
AM_CONDITIONAL(HAVE_EXTERNAL_BLAS, test x$HAVE_EXTERNAL_BLAS = xyes)
AM_CONDITIONAL(HAVE_EXTERNAL_LAPACK, test x$HAVE_EXTERNAL_LAPACK = xyes)

DO_BLASBENCH=no
if test x"$enable_debug" = x"yes" || test x"$HAVE_EXTERNAL_BLAS" = x"yes"; then
DO_BLASBENCH=yes
fi
AM_CONDITIONAL(DO_BLASBENCH, test x$DO_BLASBENCH = xyes)

DO_LAPACKBENCH=no
if test x"$enable_debug" = x"yes" || test x"$HAVE_EXTERNAL_LAPACK" = x"yes"; then
DO_LAPACKBENCH=yes
fi
AM_CONDITIONAL(DO_LAPACKBENCH, test x$DO_LAPACKBENCH = xyes)

dnl __float128 support part II.
if test x$__ICC = x"icc"; then
CFLAGS="$CFLAGS $_GCC_FLOAT128_FLAGS"
CXXFLAGS="$CXXFLAGS $_GCC_FLOAT128_FLAGS"
FCFLAGS="$FCFLAGS $_GCC_FLOAT128_FLAGS"
fi

if test x"$enable_optimization" = x"yes"; then
    CFLAGS="$CFLAGS -O2 -fomit-frame-pointer -finline-functions -funroll-all-loops"
    CXXFLAGS="$CXXFLAGS -O2 -fomit-frame-pointer -finline-functions -funroll-all-loops"
    FCFLAGS="$FCFLAGS -O2 -fomit-frame-pointer -finline-functions -funroll-all-loops"
    if test x$__ICC = x"icc"; then
        CFLAGS="$CFLAGS -xHost -O3"
        CXXFLAGS="$CXXFLAGS -xHost -O3"
        FCFLAGS="$FCFLAGS -xHost -O3"
    fi
    AC_SUBST(CFLAGS)
    AC_SUBST(CXXFLAGS)
    AC_SUBST(FCFLAGS)
fi

AC_CONFIG_FILES([
mpfrc++/Makefile
mblas/reference/Makefile
mblas/optimized/Makefile
mblas/optimized/gmp/Makefile
mblas/optimized/mpfr/Makefile
mblas/optimized/qd/Makefile
mblas/optimized/dd/Makefile
mblas/optimized/dd/cuda/Makefile
mblas/optimized/double/Makefile
mblas/optimized/__float128/Makefile
mblas/blaswrapper/Makefile
mblas/fortranwrapper/Makefile
mblas/debug/Makefile
mblas/debug/gmp/Makefile
mblas/debug/mpfr/Makefile
mblas/debug/qd/Makefile
mblas/debug/dd/Makefile
mblas/debug/double/Makefile
mblas/debug/__float128/Makefile
mblas/Makefile
mlapack/reference/Makefile
mlapack/debug/Makefile
mlapack/debug/gmp/Makefile
mlapack/debug/mpfr/Makefile
mlapack/debug/qd/Makefile
mlapack/debug/dd/Makefile
mlapack/debug/double/Makefile
mlapack/debug/__float128/Makefile
mlapack/optimized/Makefile
mlapack/optimized/gmp/Makefile
mlapack/optimized/mpfr/Makefile
mlapack/optimized/qd/Makefile
mlapack/optimized/dd/Makefile
mlapack/optimized/double/Makefile
mlapack/optimized/__float128/Makefile
mlapack/lapackwrapper/Makefile
mlapack/fortranwrapper/Makefile
mlapack/Makefile
benchmark/Makefile
benchmark/mblas/Makefile
benchmark/mlapack/Makefile
external/Makefile
external/gmp/Makefile
external/mpfr/Makefile
external/mpc/Makefile
external/qd/Makefile
external/lapack/Makefile
examples/Makefile
examples/mblas/Makefile
examples/mlapack/Makefile
])

AC_OUTPUT

