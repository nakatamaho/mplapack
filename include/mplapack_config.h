/* include/mplapack_config.h.  Generated from mplapack_config.h.in by configure.  */
/* include/mplapack_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "mplapack"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "mplapack"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "mplapack 2.0.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "mplapack"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.0.1"

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "2.0.1"

/* Define if native _Float64x is available. */
/* #undef ___MPLAPACK_HAVE_NATIVE__FLOAT64X___ */

/* Define if long double is binary128 and _Float128 is unavailable. */
/* #undef ___MPLAPACK_LONGDOUBLE_IS_BINARY128___ */

/* Define if _Float64x is emulated using __float80. */
/* #undef ___MPLAPACK_USE___FLOAT80___ */

/* Define if _Float128 uses __float128 + libquadmath. */
/* #undef ___MPLAPACK_WANT_LIBQUADMATH___ */

/* Define if _Float128 is available and equals long double. */
/* #undef ___MPLAPACK__FLOAT128_IS_LONGDOUBLE___ */

/* Define if _Float128 is available (and not long double). */
#define ___MPLAPACK__FLOAT128_ONLY___ 1

/* Define if _Float64x is compatible with long double. */
#define ___MPLAPACK__FLOAT64X_IS_LONGDOUBLE___ 1

[
/* ===== MPLAPACK public types (added by configure.ac; do not edit generated header) ===== */
#ifndef _MPLAPACK_PUBLIC_TYPES_DEFINED_
#define _MPLAPACK_PUBLIC_TYPES_DEFINED_

#ifdef __cplusplus
#include <complex>
#endif
#include <inttypes.h>
#include <stdlib.h>

/* Prefer 64-bit integer by default (same behavior as the historical mplapack_config.h). */
#ifndef USE64BITINT
#define USE64BITINT 1
#endif

#ifdef USE64BITINT
#if defined _WIN32
typedef long int mplapackint;
#elif defined __APPLE__
typedef long mplapackint;
#else
typedef int64_t mplapackint;
#endif
#endif

typedef mplapackint mplapacklogical;

#ifdef __cplusplus
typedef mplapacklogical (*LFP)(...);
#else
typedef mplapacklogical (*LFP)();
#endif

#endif /* _MPLAPACK_PUBLIC_TYPES_DEFINED_ */
]
