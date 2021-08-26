dnl
dnl Configure for MPLAPACK
dnl

m4_define([MPLAPACK_VER_MAJOR], [1])
m4_define([MPLAPACK_VER_MINOR], [0])
m4_define([MPLAPACK_VER_PATCH], [0])
m4_define([MPLAPACK_VERSION], [MPLAPACK_VER_MAJOR.MPLAPACK_VER_MINOR.MPLAPACK_VER_PATCH])
AC_INIT(mplapack, [MPLAPACK_VERSION])

AC_SUBST([PACKAGE_DESCRIPTION], "A multiple precision version of BLAS and LAPACK")
AC_CONFIG_HEADER(include/config.h) 

LAPACK_VER_MAJOR="3"
LAPACK_VER_MINOR="9"
LAPACK_VER_PATCH="1"
LAPACK_VERSION=${LAPACK_VER_MAJOR}.${LAPACK_VER_MINOR}.${LAPACK_VER_PATCH}
AC_SUBST(LAPACK_VERSION)

sed -e "s/%%mplapack_ver_major%%/MPLAPACK_VER_MAJOR/g" \
    -e "s/%%mplapack_ver_minor%%/MPLAPACK_VER_MINOR/g" \
    -e "s/%%mplapack_ver_patch%%/MPLAPACK_VER_PATCH/g" \
    -e "s/%%lapack_ver_major%%/$LAPACK_VER_MAJOR/g" \
    -e "s/%%lapack_ver_minor%%/$LAPACK_VER_MINOR/g" \
    -e "s/%%lapack_ver_patch%%/$LAPACK_VER_PATCH/g" \
     mplapack/reference/iMlaver.cpp.in > mplapack/reference/iMlaver.cpp

AC_CONFIG_FILES([ Makefile mplapack.pc ])
AC_CONFIG_FILES([packaging/PKGBUILD packaging/mplapack.spec])

AC_CANONICAL_TARGET

AC_PROG_CC
AC_PROG_CXX
IS_INTELCC=0
AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([[
    #if !defined __INTEL_COMPILER
    #error
    #endif
    #include <stdio.h>
    int main() {
        return 0;    
    }
   ]]) 
], [is_intelcc=1])
AM_CONDITIONAL(IS_INTELCC, test x$is_intelcc = x1)

AC_PROG_LIBTOOL
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_MACRO_DIR([m4])

dnl OpenMP Compiler flags.
AC_LANG([C++])
AC_OPENMP

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
dnl all mp versions of mplapack and mpblas, etc.
AC_ARG_ENABLE(test,
[  --enable-test    Enable test])
AM_CONDITIONAL(ENABLE_TEST, test x$enable_test = xyes)

dnl disable optimization and add more info to gdb.
AC_ARG_ENABLE(debug,
[  --enable-debug    Enable debug])
AM_CONDITIONAL(ENABLE_DEBUG, test x$enable_debug = xyes)

AC_ARG_ENABLE(optimization,
AS_HELP_STRING([--enable-optimization], [Enable compiler level optimization]),
[case "${enableval}" in
    yes) enable_optimization=yes ;;
    no)  enable_optimization=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-optimization) ;;
        esac],[enable_optimization=yes])

if test x"$enable_optimization" = x"yes" && test x"$enable_debug" = x"yes" ; then
        AC_MSG_ERROR([you cannot set both enable-optimization=yes and enable-debug=yes])
fi

AC_ARG_ENABLE(cuda,
[  --enable-cuda    Enable compiler level cuda])
AM_CONDITIONAL(ENABLE_CUDA, test x$enable_cuda = xyes)

AC_ARG_WITH(cudatoolkithome,
AS_HELP_STRING([--with-cudatoolkithome], [Specify path to CUDA Toolkit (default: /usr/local/cuda)]),
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

dnl disable benchmark
AC_ARG_ENABLE(benchmark,
AS_HELP_STRING([--disable-benchmark], [Disable benchmark]))
AM_CONDITIONAL(DISABLE_BENCHMARK, test x$enable_benchmark = xno)

dnl disable examples
AC_ARG_ENABLE(examples, 
AS_HELP_STRING([--disable-examples], [Disable examples]))
AM_CONDITIONAL(DISABLE_EXAMPLES, test x$enable_examples = xno)

dnl disable reference
AC_ARG_ENABLE(reference,
AS_HELP_STRING([--disable-reference], [Build without reference libraries]))
AM_CONDITIONAL(DISABLE_REFERENCE, test x$enable_reference = xno)

dnl enable GMP version of MPLAPACK.
AC_ARG_ENABLE(gmp,
AS_HELP_STRING([--enable-gmp], [Enable GMP version]))
AM_CONDITIONAL(ENABLE_GMP, test x$enable_gmp = xyes)

dnl enable MPFR version of MPLAPACK.
AC_ARG_ENABLE(mpfr,
AS_HELP_STRING([--enable-mpfr], [Enable MPFR version]))

dnl We use MPFR version as reference
dnl so if test is turned on, 
if test x$enable_test = x"yes"; then
    enable_mpfr=yes
fi
AM_CONDITIONAL(ENABLE_MPFR, test x$enable_mpfr = xyes)

dnl enable QD version of MPLAPACK.
AC_ARG_ENABLE(qd,
AS_HELP_STRING([--enable-qd], [Enable QD version]))
AM_CONDITIONAL(ENABLE_QD, test x$enable_qd = xyes)

AC_ARG_WITH(qdlib-flags,
AS_HELP_STRING([--with-qdlib-flags], [passing flag to internal qd library])
,,with_qdlib_flags="")
QDLIB_FLAGS=$with_qdlib_flags;
AC_SUBST(QDLIB_FLAGS)

dnl enable DD version of MPLAPACK.
AC_ARG_ENABLE(dd,
AS_HELP_STRING([--enable-dd], [Enable DD version]),
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
AS_HELP_STRING([--enable-double], [Enable double version]),
[case "${enableval}" in
    yes) enable_double=yes ;;
    no)  enable_double=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-double) ;;
        esac],[enable_double=no])
AM_CONDITIONAL(ENABLE_DOUBLE, test x$enable_double = xyes)

AC_ARG_ENABLE(_Float128,
AS_HELP_STRING([--enable-_Float128], [Enable _Float128 version (TS 18661-3 like gcc extension) ]),
[case "${enableval}" in
    yes) enable__Float128=yes ;;
    no)  enable__Float128=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-_Float128) ;;
        esac],[enable__Float128=no])

SAVE_CXXFLAGS="$CXXFLAGS"
AM_CONDITIONAL(WANT_QUADMATH, test x$enable__Float128 = xno)
if test x$enable__Float128 = x"yes"; then
# For unknown reason, the f128 suffix doesn't work in g++-9, but it does in gcc-9.
for gcc works.
    AC_MSG_CHECKING([whether _Float128 is supported])
    AC_LINK_IFELSE([
       AC_LANG_SOURCE([[
       #include <stdio.h>
       #include <math.h>
       #include <stdlib.h>
       int main() {
           _Float128 a, b, c; b = 1.0q; c=3.0q; a = b / c;
           char buf[1024];
           strfromf128(buf, sizeof(buf), "%.40g", a);
           printf("fp128 %s\n",buf);
       }
      ]])
    ], [have__Float128=yes want_libquadmath=no], [have__Float128=no])

    if test x"$have__Float128" = x"yes"; then    
        AC_COMPILE_IFELSE([
        AC_LANG_SOURCE([[
           #include <stdio.h>
           #include <math.h>
           #include <stdlib.h>
           int main() {
               long double *x;
               _Float128 *y = x;
           }
        ]]) 
        ], [_Float128islongdouble=yes], [_Float128islongdouble=no])

        if test "$_Float128islongdouble" = "yes"; then 
             AC_MSG_RESULT([_Float128 is long double])
        else
             AC_MSG_RESULT([_Float128 is not long double])
        fi
    fi

    if test x"$have__Float128" = x"no"; then
        AC_COMPILE_IFELSE([
        AC_LANG_SOURCE([[
           #include <stdio.h>
           #include <limits>
           int main(){
           #if __LDBL_MANT_DIG__  >= 113
              printf("%d\n",__LDBL_MANT_DIG__);
           #else
           #error
           #endif
        }
        ]]) 
        ], [longdoubleisbinary128=yes], [longdoubleisbinary128=no])

        if test "$longdoubleisbinary128" = "yes"; then 
             AC_MSG_RESULT([long double is binary128])
        else
             AC_MSG_RESULT([long double is not binary128])
        fi

        save_LIBS="$LIBS"
        LIBS="-lquadmath $LIBS"
        AC_LINK_IFELSE([
            AC_LANG_SOURCE([[
                #include <quadmath.h>
                #include <stdio.h>
                int main(){
                    char buf[1024];
                    __float128 a, b, c; b = 1.0q; c= 3.0q; a = b / c;
                    quadmath_snprintf(buf,128,"%.40Qf\n",a);
                    printf("%s",buf);
                }
           ]])
        ], [have___float128=yes want_libquadmath=yes], [have___float128=no])
        LIBS="$save_LIBS"
    fi

    if test "$have___float128" = "yes" && test x"$want_libquadmath" = x"yes" ; then
        AC_MSG_RESULT([_Float128 is supported via __float128 and libquadmath])
        sed 's/%%MPLAPACK_BINARY128%%/#define ___MPLAPACK_WANT_LIBQUADMATH___/g' include/mplapack_utils__Float128.h.in > include/mplapack_utils__Float128.h
    elif test "$have__Float128" = "yes" && test x"$want_libquadmath" = x"no" && test x"$_Float128islongdouble" = x"yes"; then
        AC_MSG_RESULT([_Float128 is supported via native libc and _Float128 and long double are the same])
        sed 's/%%MPLAPACK_BINARY128%%/#define ___MPLAPACK__FLOAT128_IS_LONGDOUBLE___/g' include/mplapack_utils__Float128.h.in > include/mplapack_utils__Float128.h
    elif test "$have__Float128" = "yes" && test x"$want_libquadmath" = x"no" && test x"$_Float128islongdouble" != x"yes"; then
        AC_MSG_RESULT([_Float128 is supported via native libc and only binary128 in this system])
        sed 's/%%MPLAPACK_BINARY128%%/#define ___MPLAPACK__FLOAT128_ONLY___/g' include/mplapack_utils__Float128.h.in > include/mplapack_utils__Float128.h
    elif test "$have__Float128" = "no" && test x"$longdoubleisbinary128" = x"yes"; then
        AC_MSG_RESULT([long double is binary128, but neither _Float128 nor __float128 are supported])
        sed 's/%%MPLAPACK_BINARY128%%/#define ___MPLAPACK_LONGDOUBLE_IS_BINARY128___/g' include/mplapack_utils__Float128.h.in > include/mplapack_utils__Float128.h
    else 
        AC_MSG_ERROR([binary128 is not supported on this system])
    fi
fi
AM_CONDITIONAL(ENABLE__FLOAT128, test x$enable__Float128 = xyes)
AM_CONDITIONAL(WANT_QUADMATH, test x$want_libquadmath = xyes)
AM_CONDITIONAL(IS__FLOAT128_LONGDOUBLE, test x$_Float128islongdouble = xyes)
AM_CONDITIONAL(IS_LONGDOUBLE_BINARY128, test x$longdoubleisbinary128 = xyes)

CXXFLAGS="$SAVE_CXXFLAGS"

AC_ARG_ENABLE(_Float64x,
AS_HELP_STRING([--enable-_Float64x], [Enable _Float64x (aka __float80 or long double) version]),
[case "${enableval}" in
    yes) enable__Float64x=yes ;;
    no)  enable__Float64x=no  ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-_Float64x) ;;
        esac],[enable__Float64x=no])

SAVE_CXXFLAGS="$CXXFLAGS"
if test x$enable__Float64x = x"yes"; then
    AC_MSG_CHECKING([whether _Float64x is supported])
    AC_LINK_IFELSE([
       AC_LANG_SOURCE([[
       #include <stdio.h>
       #include <math.h>
       #include <stdlib.h>
       int main() {
           _Float64x a, b, c; b = 1.0w; c=3.0w; a = b / c;
           printf("_Float64x %Lf\n", a);
       }
      ]])
    ], [have__Float64x=yes want___float80=no], [_have__Float64x=no])

    if test "$_have__Float64x" = "no"; then
        AC_LINK_IFELSE([
            AC_LANG_SOURCE([[
                #include <stdio.h>
                #include <math.h>
                #include <stdlib.h>
                int main(){
                    char buf[1024];
                    __float80 a, b, c; b = 1.0w; c= 3.0w; a = b / c;
                    printf("__float80 %Lf",a);
                }
           ]])
        ], [have__Float64x=yes want___float80=yes], [have__Float64x=no])

        if test "$have__Float64x" = "no"; then
            AC_MSG_ERROR([_Float64x is not supported on this platform])
        fi
        LDFLAGS="$save_LDFLAGS"
    fi

    if test x"$have__Float64x" = x"yes" && test x"$want___float80" = x"no" ; then    
        AC_MSG_RESULT([_Float64x is supported natively])
        AC_COMPILE_IFELSE([
        AC_LANG_SOURCE([[
           #include <stdio.h>
           #include <math.h>
           #include <stdlib.h>
           int main() {
               long double *x;
               _Float64x *y = x;
           }
        ]]) 
        ], [float64xislongdouble=yes], [float64islongdouble=no])
    fi
    if test x"$have__Float64x" = x"yes" && test x"$want___float80" = x"yes" ; then    
        AC_MSG_RESULT([_Float64x is supported via __float80])
        AC_COMPILE_IFELSE([
        AC_LANG_SOURCE([[
           #include <stdio.h>
           #include <math.h>
           #include <stdlib.h>
           int main() {
               long double *x;
               __float80 *y = x;
           }
        ]]) 
        ], [float64xislongdouble=yes], [float64islongdouble=no])
    fi
    if test "$float64xislongdouble" = "yes"; then 
         AC_MSG_RESULT([_Float64x is long double])
    else
         AC_MSG_RESULT([_Float64x is not long double])
    fi
    if test "$want___float80" = "yes" && test "$float64xislongdouble" = "yes" ; then
        sed 's/%%MPLAPACK__FLOAT64X___FLOAT80%%/typedef __float80 _Float64x;\n#define ___MPLAPACK__FLOAT64X_IS_LONGDOUBLE___/g' include/mplapack_utils__Float64x.h.in > include/mplapack_utils__Float64x.h
    fi
    if test "$want___float80" = "yes" && test "$float64xislongdouble" = "no" ; then
        sed 's/%%MPLAPACK__FLOAT64X___FLOAT80%%/typedef __float80 _Float64x;/g' include/mplapack_utils__Float64x.h.in > include/mplapack_utils__Float64x.h
    fi
    if test "$want___float80" = "no" && test "$float64xislongdouble" = "yes" ; then
        sed 's/%%MPLAPACK__FLOAT64X___FLOAT80%%/#define ___MPLAPACK__FLOAT64X_IS_LONGDOUBLE___/g' include/mplapack_utils__Float64x.h.in > include/mplapack_utils__Float64x.h
    fi
    if test "$want___float80" = "no" && test "$float64xislongdouble" = "no" ; then
        sed 's/%%MPLAPACK__FLOAT64X___FLOAT80%%//g' include/mplapack_utils__Float64x.h.in > include/mplapack_utils__Float64x.h
    fi
fi
CXXFLAGS="$SAVE_CXXFLAGS"

AM_CONDITIONAL(ENABLE__FLOAT64X, test x$enable__Float64x = xyes)
AM_CONDITIONAL(IS__FLOAT64X_LONGDOUBLE, test x$float64xislongdouble = xyes)
AM_CONDITIONAL(IS___FLOAT80, test x$want___float80 = xyes)

highest_abs_top_srcdir=`cd "$srcdir" && pwd` 

dnl GMP library start
AC_ARG_WITH(system-gmp,
AS_HELP_STRING([--with-system-gmp=PREFIX],
	    [ ="yes" if it is at the default location (/usr/local/lib), or alternative prefix path to gmp library]
	  ),
	  GMP_PATHSET=1,
	  GMP_PATHSET=0;
	  with_system_gmp="no";
	  )
GMP_INCLUDEDIR=$highest_abs_top_srcdir/external/i/GMP/include/
GMP_LIBDIR=$highest_abs_top_srcdir/external/i/GMP/lib/
if test  $GMP_PATHSET = 1 ; then
	if test $with_system_gmp != "yes"; then
		GMP_INCLUDEDIR="$with_system_gmp/include"
		GMP_LIBDIR="$with_system_gmp/lib"
	fi
fi

if test x$enable_gmp = x"yes" -o x$enable_mpfr = x"yes" -o x$enable_test = x"yes"; then
    if test x$with_system_gmp != x"no"; then
        HAVE_GMP=""
        save_CXXFLAGS=$CXXFLAGS
        CXXFLAGS="$CXXFLAGS -I$GMP_INCLUDEDIR"
        AC_CHECK_HEADER(gmp.h, HAVE_GMP_HEADER=yes)
        CXXFLAGS=$save_CXXFLAGS
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
        AC_MSG_RESULT([Using GMP library comes with MPLAPACK])
        GMP_INCLUDEDIR=$highest_abs_top_srcdir/external/i/GMP/include/
        GMP_LIBDIR=$highest_abs_top_srcdir/external/i/GMP/lib/
    fi
fi
AM_CONDITIONAL(SYSTEM_GMP, test x$SYSTEM_GMP = xyes)
AC_SUBST(GMP_INCLUDEDIR, [$GMP_INCLUDEDIR])
AC_SUBST(GMP_LIBDIR, [$GMP_LIBDIR])
dnl GMP library end

dnl MPFR library start
AC_ARG_WITH(system-mpfr,
AS_HELP_STRING([--with-system-mpfr=PREFIX],
	    [ ="yes" if it is at the default location (/usr/local/lib), or alternative prefix path to mpfr library]
	  ),
	  MPFR_PATHSET=1,
	  MPFR_PATHSET=0;
	  with_system_mpfr="no";
	  )
MPFR_INCLUDEDIR=$highest_abs_top_srcdir/external/i/MPFR/include/
MPFR_LIBDIR=$highest_abs_top_srcdir/external/i/MPFR/lib/
if test  $MPFR_PATHSET = 1 ; then
	if test $with_system_mpfr != "yes"; then
		MPFR_INCLUDEDIR="$with_system_mpfr/include"
		MPFR_LIBDIR="$with_system_mpfr/lib"
	fi
fi
if test x$enable_mpfr = x"yes" -o x$enable_test = x"yes"; then
    if test x$with_system_mpfr != x"no"; then
        HAVE_MPFR=""
        save_CXXFLAGS=$CXXFLAGS
        CXXFLAGS="$CXXFLAGS -I$MPFR_INCLUDEDIR -I$GMP_INCLUDEDIR"
        AC_CHECK_HEADER(mpfr.h, HAVE_MPFR_HEADER=yes)
        if test x$HAVE_MPFR_HEADER = xyes; then
            save_LDFLAGS="$LDFLAGS"
            LDFLAGS="$LDFLAGS -L$MPFR_LIBDIR -L$GMP_LIBDIR"
            AC_CHECK_LIB(mpfr, mpfr_init, [MPFR_LIBS="-lmpfr -lgmp"; HAVE_MPFR=yes])
            LDFLAGS=$save_LDFLAGS
        fi
        if test x$HAVE_MPFR != xyes; then
            AC_MSG_RESULT([No MPFR library found; using internal one instead])
        else
            AC_MSG_RESULT([MPFR found: using MPFR comes with the system])
            SYSTEM_MPFR=yes
        fi
        CXXFLAGS=$save_CXXFLAGS
    else
        AC_MSG_RESULT([Using MPFR comes with MPLAPACK])
        MPFR_INCLUDEDIR=$highest_abs_top_srcdir/external/i/MPFR/include/
        MPFR_LIBDIR=$highest_abs_top_srcdir/external/i/MPFR/lib/
    fi
fi
AM_CONDITIONAL(SYSTEM_MPFR, test x$SYSTEM_MPFR = xyes)
AC_SUBST(MPFR_INCLUDEDIR, [$MPFR_INCLUDEDIR])
AC_SUBST(MPFR_LIBDIR, [$MPFR_LIBDIR])
dnl MPFR library end

dnl MPC library start
AC_ARG_WITH(system-mpc,
AS_HELP_STRING([--with-system-mpc=PREFIX],
	    [ ="yes" if it is at the default location (/usr/local/lib), or alternative prefix path to mpc library]
	  ),
	  MPC_PATHSET=1,
	  MPC_PATHSET=0;
	  with_system_mpc="no";
	  )
MPC_INCLUDEDIR=$highest_abs_top_srcdir/external/i/MPC/include/
MPC_LIBDIR=$highest_abs_top_srcdir/external/i/MPC/lib/
if test  $MPC_PATHSET = 1 ; then
	if test $with_system_mpc != "yes"; then
		MPC_INCLUDEDIR="$with_system_mpc/include"
		MPC_LIBDIR="$with_system_mpc/lib"
	fi
fi
if test x$enable_mpfr = x"yes" -o x$enable_test = x"yes"; then
    if test x$with_system_mpc != x"no"; then
        HAVE_MPC=""
        save_CXXFLAGS=$CXXFLAGS
        CXXFLAGS="$CXXFLAGS -I$MPC_INCLUDEDIR -I$MPFR_INCLUDEDIR -I$GMP_INCLUDEDIR"
        AC_CHECK_HEADER(mpc.h, HAVE_MPC_HEADER=yes)
        if test x$HAVE_MPC_HEADER = xyes; then
            save_LDFLAGS="$LDFLAGS"
            LDFLAGS="$LDFLAGS -L$MPC_LIBDIR -L$MPFR_LIBDIR -L$GMP_LIBDIR"
            AC_CHECK_LIB(mpc, mpc_init2, [MPC_LIBS="-lmpc -lmpfr -lgmp"; HAVE_MPC=yes])
            LDFLAGS=$save_LDFLAGS
        fi
        if test x$HAVE_MPC != xyes; then
            AC_MSG_RESULT([No MPC library found; using internal one instead])
        else
            AC_MSG_RESULT([MPC library found: using MPC library comes with the system])
            SYSTEM_MPC=yes
        fi
        CXXFLAGS=$save_CXXFLAGS
    else
	AC_MSG_RESULT([Using MPC library comes with MPLAPACK])
        MPC_INCLUDEDIR=$highest_abs_top_srcdir/external/i/MPC/include/
        MPC_LIBDIR=$highest_abs_top_srcdir/external/i/MPC/lib/
    fi
fi
AM_CONDITIONAL(SYSTEM_MPC, test x$SYSTEM_MPC = xyes)
AC_SUBST(MPC_INCLUDEDIR, [$MPC_INCLUDEDIR])
AC_SUBST(MPC_LIBDIR, [$MPC_LIBDIR])
dnl MPC library end

dnl QD library start
AC_ARG_WITH(system-qd,
AS_HELP_STRING([--with-system-qd=PREFIX],
	    [ ="yes" if it is at the default location (/usr/local/lib), or alternative prefix path to qd library]
	  ),
	  QD_PATHSET=1,
	  QD_PATHSET=0;
	  with_system_qd="no";
	  )
QD_INCLUDEDIR=$highest_abs_top_srcdir/external/i/QD/include/
QD_LIBDIR=$highest_abs_top_srcdir/external/i/QD/lib/
if test  $QD_PATHSET = 1 ; then
	if test $with_system_qd != "yes"; then
		QD_INCLUDEDIR="$with_system_qd/include"
		QD_LIBDIR="$with_system_qd/lib"
	fi
fi
if test x$enable_qd = x"yes" -o x$enable_dd = x"yes" -o x$enable_test = x"yes"; then
    if test x$with_system_qd != x"no"; then
        HAVE_QD=""
        save_CXXFLAGS=$CXXFLAGS
        CXXFLAGS="$CXXFLAGS -I$QD_INCLUDEDIR"
        AC_CHECK_HEADER(qd/qd_real.h, HAVE_QD_HEADER=yes)
        CXXFLAGS=$save_CXXFLAGS
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
        fi
    else
        AC_MSG_RESULT([Using QD library comes with MPLAPACK])
        QD_INCLUDEDIR=$highest_abs_top_srcdir/external/i/QD/include/
        QD_LIBDIR=$highest_abs_top_srcdir/external/i/QD/lib/
    fi   
fi
AM_CONDITIONAL(SYSTEM_QD, test x$SYSTEM_QD = xyes)
AC_SUBST(QD_INCLUDEDIR, [$QD_INCLUDEDIR])
AC_SUBST(QD_LIBDIR, [$QD_LIBDIR])
dnl QD library end

need_internal_gmplib=no
if test x$enable_gmp = x"yes" || test x$enable_mpfr = x"yes" && test x$SYSTEM_GMP != x"yes"; then
    need_internal_gmplib=yes
    AC_SUBST(GMP_DIR, [gmp])
fi
#AM_CONDITIONAL(NEED_INTERNAL_GMPLIB, test x$need_internal_gmplib = xyes)

need_internal_qdlib=no
if test x$enable_qd = x"yes" || test x$enable_dd = x"yes"; then
    if test x$SYSTEM_QD != x"yes"; then
       need_internal_qdlib=yes
       AC_SUBST(QD_DIR, [qd])
    fi
fi
#AM_CONDITIONAL(NEED_INTERNAL_QDLIB, test x$need_internal_qdlib = xyes)

need_internal_mpfrlib=no
if test x$enable_mpfr = x"yes" && test x$SYSTEM_MPFR != x"yes"; then
    need_internal_mpfrlib=yes
    AC_SUBST(MPFR_DIR, [mpfr])
fi
#AM_CONDITIONAL(NEED_INTERNAL_MPFRLIB, test x$need_internal_mpfrlib = xyes)

need_internal_mpclib=no
if test x$enable_mpfr = x"yes" && test x$SYSTEM_MPC != x"yes"; then
    need_internal_mpclib=yes
    AC_SUBST(MPC_DIR, [mpc])
fi
#AM_CONDITIONAL(NEED_INTERNAL_MPCLIB, test x$need_internal_mpclib = xyes)

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

if test x$enable_qd = x"yes"; then
AC_MSG_RESULT([Enable QD version])
else
AC_MSG_RESULT([Disable QD version])
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
if test x$enable__Float64x = x"yes"; then
    if test "$want___float80" = "yes"; then
        AC_MSG_RESULT([Enable _Float64x version via __float80])
    else
        AC_MSG_RESULT([Enable native _Float64x version ])
    fi
else
    AC_MSG_RESULT([Disable _Float64x version])
fi
if test x$enable__Float64x = x"yes"; then
    if test "$float64xislongdouble" = "yes"; then
        AC_MSG_RESULT([_Float64x is long double])
    else
        AC_MSG_RESULT([_Float64x is not long double])
    fi
fi

if test x$enable__Float128 = x"yes" && test x"$want_libquadmath" = x"yes"; then
        AC_MSG_RESULT([_Float128 is supported via __float128 and libquadmath])
elif test x$enable__Float128 = x"yes" && test x"$want_libquadmath" = x"no" && test x"$_Float128islongdouble" = x"yes"; then
        AC_MSG_RESULT([_Float128 is supported via native libc and _Float128 and long double are the same])
elif test x$enable__Float128 = x"yes" && test x"$want_libquadmath" = x"no" && test x"$_Float128islongdouble" != x"yes"; then
        AC_MSG_RESULT([_Float128 is supported via native libc and only binary128 in this system])
elif test x$enable__Float128 = x"yes" && test x"$longdoubleisbinary128" = x"yes"; then
        AC_MSG_RESULT([long double is binary128, but neither _Float128 nor __float128 are supported])
else
        AC_MSG_RESULT([disable _Float128])
fi

if test x$enable_test = x"yes"; then
AC_MSG_RESULT([Enable test])
AC_SUBST(LAPACK_DIR, [lapack])
LAPACK_INCLUDEDIR=$highest_abs_top_srcdir/external/i/LAPACK/include/
LAPACK_LIBDIR=$highest_abs_top_srcdir/external/i/LAPACK/lib/
AC_SUBST(LAPACK_INCLUDEDIR, [$LAPACK_INCLUDEDIR])
AC_SUBST(LAPACK_LIBDIR, [$LAPACK_LIBDIR])
else
LAPACK_INCLUDEDIR=.
LAPACK_LIBDIR=.
AC_SUBST(LAPACK_INCLUDEDIR, [$LAPACK_INCLUDEDIR])
AC_SUBST(LAPACK_LIBDIR, [$LAPACK_LIBDIR])
AC_MSG_RESULT([Disable test])
fi

dnl print status...end

AC_ARG_WITH(external-blas,
AS_HELP_STRING([--with-external-blas], [Specify location of external BLAS e.g. "-L/usr/local/lib -lblas". It is used for benchmark only.]))

AC_ARG_WITH(external-lapack,
AS_HELP_STRING([--with-external-lapack], [Specify how to link against external LAPACK e.g. "-L/usr/local/lib -llapack". It is used for benchmark only.]))

dnl First we check Fotran 90 compiler stuffs.

if test x"$enable_test" = x"yes"; then
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
if test x"$enable_test" = x"yes" || test x"$HAVE_EXTERNAL_BLAS" = x"yes"; then
DO_BLASBENCH=yes
fi
AM_CONDITIONAL(DO_BLASBENCH, test x$DO_BLASBENCH = xyes)

DO_LAPACKBENCH=no
if test x"$enable_test" = x"yes" || test x"$HAVE_EXTERNAL_LAPACK" = x"yes"; then
DO_LAPACKBENCH=yes
fi
AM_CONDITIONAL(DO_LAPACKBENCH, test x$DO_LAPACKBENCH = xyes)

if test x"$enable_optimization" = x"no" && test x"$enable_debug" = x"yes" ; then
    CFLAGS="-O0 -g"
    CXXFLAGS="-O0 -g"
    FCFLAGS="-O0 -g"
fi

if test x"$enable_optimization" = x"yes" && test x"$enable_debug" = x"no" ; then    
    case $target in
       powerpc*-*-*)
           CFLAGS="$CFLAGS -mcpu=native -O2 -finline-functions -funroll-all-loops"
           CXXFLAGS="$CXXFLAGS -mcpu=native -O2 -finline-functions -funroll-all-loops"
           FCFLAGS="$FCFLAGS -mcpu=native -O2 -finline-functions -funroll-all-loops"
       ;;
       * ) 
           CFLAGS="$CFLAGS -march=native -g -O2 -finline-functions -funroll-all-loops"
           CXXFLAGS="$CXXFLAGS -march=native -g -O2 -finline-functions -funroll-all-loops"
           FCFLAGS="$FCFLAGS -march=native -g -O2 -finline-functions -funroll-all-loops"
       ;;
    esac
    if test x"is_intelcc" = x"1"; then
        CXXFLAGS="$CXXFLAGS -xHost -O3"
        CFLAGS="$CFLAGS -xHost -O3"
        FCFLAGS="$FCFLAGS -xHost -O3"
    fi
fi

AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(FCFLAGS)

############ select created libraries for the pkg-config ###############

PACKAGE_CREATED_LIBS=""
if test x$enable_dd = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack_dd -lmpblas_dd"
fi
if test x$enable_qd = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack_qd -lmpblas_qd"
fi
if test x$enable_gmp = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack_gmp -lmpblas_gmp"
fi
if test x$enable_mpfr = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack_mpfr -lmpblas_mpfr"
fi
if test x$enable_mpc = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack_mpc -lmpblas_mpc"
fi
if test x$enable__Float128 = x"yes"; then
    if test "$want_libquadmath" = "yes"; then
        PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack__Float128 -lmpblas__Float128 -lquadmath"
    else
        PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack__Float128 -lmpblas__Float128"
    fi
fi
if test x$enable_double = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack_double -lmpblas_double"
fi
if test x$enable__Float64x = x"yes"; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmplapack__Float64x -lmpblas__Float64x"
fi

# add external libraries only if they will be installed with the package
if test x$need_internal_gmplib = xyes; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lgmpxx -lgmp"
fi
if test x$need_internal_qdlib = xyes; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lqd"
fi
if test x$need_internal_mpfrlib = xyes; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmpfr"
fi
if test x$need_internal_mpclib = xyes; then
  PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" -lmpc"
fi

PACKAGE_CREATED_LIBS=$PACKAGE_CREATED_LIBS" ${OPENMP_CXXFLAGS}"

AC_SUBST([PACKAGE_CREATED_LIBS])

########################################################################


AC_CONFIG_FILES([
mpblas/optimized/Makefile
mpblas/optimized/gmp/Makefile
mpblas/optimized/mpfr/Makefile
mpblas/optimized/qd/Makefile
mpblas/optimized/dd/Makefile
mpblas/optimized/dd/cuda/Makefile
mpblas/optimized/double/Makefile
mpblas/optimized/_Float64x/Makefile
mpblas/optimized/_Float128/Makefile
mpblas/test/Makefile
mpblas/test/gmp/Makefile
mpblas/test/mpfr/Makefile
mpblas/test/qd/Makefile
mpblas/test/dd/Makefile
mpblas/test/double/Makefile
mpblas/test/_Float64x/Makefile
mpblas/test/_Float128/Makefile
mpblas/Makefile
mplapack/test/Makefile
mplapack/test/compare/Makefile
mplapack/test/compare/gmp/Makefile
mplapack/test/compare/mpfr/Makefile
mplapack/test/compare/qd/Makefile
mplapack/test/compare/dd/Makefile
mplapack/test/compare/double/Makefile
mplapack/test/compare/_Float64x/Makefile
mplapack/test/compare/_Float128/Makefile
mplapack/test/matgen/Makefile
mplapack/test/lin/Makefile
mplapack/test/lin/mpfr/Makefile
mplapack/test/lin/gmp/Makefile
mplapack/test/lin/qd/Makefile
mplapack/test/lin/dd/Makefile
mplapack/test/lin/double/Makefile
mplapack/test/lin/_Float64x/Makefile
mplapack/test/lin/_Float128/Makefile
mplapack/test/eig/Makefile
mplapack/test/eig/mpfr/Makefile
mplapack/test/eig/gmp/Makefile
mplapack/test/eig/qd/Makefile
mplapack/test/eig/dd/Makefile
mplapack/test/eig/double/Makefile
mplapack/test/eig/_Float64x/Makefile
mplapack/test/eig/_Float128/Makefile
mplapack/Makefile
external/Makefile
external/gmp/Makefile
external/mpfr/Makefile
external/mpc/Makefile
external/qd/Makefile
external/lapack/Makefile
])

if test x$enable_benchmark != xno; then
AC_CONFIG_FILES([
benchmark/Makefile
benchmark/mpblas/Makefile
benchmark/mplapack/Makefile
])
fi
if test x$enable_examples != xno; then
AC_CONFIG_FILES([
examples/Makefile
examples/mpblas/Makefile
examples/mplapack/Makefile
examples/mplapack/04_NonsymmetricEigenproblems/Makefile
examples/mplapack/05_SingularValueDecomposition/Makefile
examples/mplapack/99_Miscellaneous/Makefile
])
fi
if test x$enable_reference != xno; then
AC_CONFIG_FILES([
mpblas/reference/Makefile
mplapack/reference/Makefile
])
fi

AC_OUTPUT

