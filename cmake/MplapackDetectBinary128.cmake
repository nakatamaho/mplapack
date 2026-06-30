# MplapackDetectBinary128.cmake
#
# Reproduces the binary128 capability probing from configure.ac. After
# inclusion (with MPLAPACK_ENABLE_BINARY128 ON) the following cache-independent
# variables are set in the caller's scope, ready to substitute into
# mplapack_config.h.in:
#
#   MPLAPACK_HAVE__FLOAT128            MPLAPACK_HAVE___FLOAT128
#   MPLAPACK_LDBL_IS_BINARY128         MPLAPACK__FLOAT128_IS_LDBL
#   MPLAPACK_FLOAT128_SAME_TYPES       MPLAPACK_FLOAT128_CONVERT_GNU_TO_STD
#   MPLAPACK_FLOAT128_CONVERT_STD_TO_GNU
#   MPLAPACK_BINARY128_MODE  (0 DISABLED / 1 LDBL / 2 FLOAT128 / 3 QUADMATH)
#   MPLAPACK_BINARY128_IO    (0 NONE / 1 SNPRINTF_LDBL / 2 STRFROMF128 / 3 QUADMATH_SNPRINTF)
#   MPLAPACK_BINARY128_MATH  (0 NONE / 1 LDBL / 2 F128 / 3 QUADMATH)
#   MPLAPACK_HAVE_STD_ABS_FLOAT128
#   MPLAPACK_HAVE_STD_MATH_FLOAT128
#   MPLAPACK_HAVE_STD_COMPLEX_FLOAT128
#   MPLAPACK_HAVE_C_COMPLEX_FLOAT128
#   MPLAPACK_BINARY128_EXTRA_LIBS / MPLAPACK_BINARY128_EXTRA_FLAGS
#
# Priority follows configure.ac: _Float128 > __float128+quadmath > long double.

include(CheckCXXSourceCompiles)

function(_mp_try var src)
  # Compile (and link) a C++ snippet; sets ${var} to 1/0 in the parent scope.
  check_cxx_source_compiles("${src}" ${var}_RESULT)
  if(${var}_RESULT)
    set(${var} 1 PARENT_SCOPE)
  else()
    set(${var} 0 PARENT_SCOPE)
  endif()
endfunction()

# Defaults (binary128 disabled).
set(MPLAPACK_HAVE__FLOAT128 0)
set(MPLAPACK_HAVE___FLOAT128 0)
set(MPLAPACK_LDBL_IS_BINARY128 0)
set(MPLAPACK__FLOAT128_IS_LDBL 0)
set(MPLAPACK_FLOAT128_SAME_TYPES 0)
set(MPLAPACK_FLOAT128_CONVERT_GNU_TO_STD 0)
set(MPLAPACK_FLOAT128_CONVERT_STD_TO_GNU 0)
set(MPLAPACK_BINARY128_MODE 0)
set(MPLAPACK_BINARY128_IO 0)
set(MPLAPACK_BINARY128_MATH 0)
set(MPLAPACK_HAVE_STD_ABS_FLOAT128 0)
set(MPLAPACK_HAVE_STD_MATH_FLOAT128 0)
set(MPLAPACK_HAVE_STD_COMPLEX_FLOAT128 0)
set(MPLAPACK_HAVE_C_COMPLEX_FLOAT128 0)
set(MPLAPACK_BINARY128_EXTRA_LIBS "")
set(MPLAPACK_BINARY128_EXTRA_FLAGS "")

if(NOT MPLAPACK_ENABLE_BINARY128)
  return()
endif()

message(STATUS "Detecting binary128 support")

# --- Phase 1: type availability -------------------------------------------
_mp_try(_mp_have_f128 "int main(){ _Float128 a,b,c; b=1.0; c=3.0; a=b/c; (void)a; return 0; }")
_mp_try(_mp_have_q128 "int main(){ __float128 a,b,c; b=1.0; c=3.0; a=b/c; (void)a; return 0; }")
_mp_try(_mp_ldbl_b128 "int main(){\n#if (__LDBL_MANT_DIG__ >= 113) && (__LDBL_MAX_EXP__ >= 16384)\nreturn 0;\n#else\n#error not binary128\n#endif\n}")

set(_mp_f128_is_ldbl 0)
if(_mp_have_f128)
  _mp_try(_mp_f128_is_ldbl "int main(){ long double *x=0; _Float128 *y=x; (void)y; return 0; }")
endif()

set(_mp_same 0)
set(_mp_gnu2std 0)
set(_mp_std2gnu 0)
if(_mp_have_f128 AND _mp_have_q128)
  _mp_try(_mp_same "#include <type_traits>\nint main(){ static_assert(std::is_same<_Float128,__float128>::value,\"\"); return 0; }")
  _mp_try(_mp_gnu2std "#include <type_traits>\nint main(){ static_assert(std::is_convertible<__float128,_Float128>::value,\"\"); return 0; }")
  _mp_try(_mp_std2gnu "#include <type_traits>\nint main(){ static_assert(std::is_convertible<_Float128,__float128>::value,\"\"); return 0; }")
endif()

# --- Phase 2: I/O functions ------------------------------------------------
set(_mp_strfromf128 0)
if(_mp_have_f128)
  _mp_try(_mp_strfromf128 "#include <stdlib.h>\nint main(){ _Float128 a=1.0; char b[128]; strfromf128(b,sizeof(b),\"%.40g\",a); return 0; }")
endif()

set(_mp_quadmath_snprintf 0)
set(_mp_quadmath_math 0)
if(_mp_have_q128)
  set(_save_libs ${CMAKE_REQUIRED_LIBRARIES})
  list(APPEND CMAKE_REQUIRED_LIBRARIES quadmath)
  _mp_try(_mp_quadmath_snprintf "#include <quadmath.h>\nint main(){ __float128 a=1.0q; char b[128]; quadmath_snprintf(b,sizeof(b),\"%.40Qf\",a); return 0; }")
  _mp_try(_mp_quadmath_math "#include <quadmath.h>\nint main(){ __float128 x=2.0q; (void)sqrtq(x);(void)expq(x);(void)logq(x);(void)sinq(x);(void)cosq(x); return 0; }")
  set(CMAKE_REQUIRED_LIBRARIES ${_save_libs})
endif()

set(_mp_snprintf_Lg 0)
if(_mp_ldbl_b128)
  _mp_try(_mp_snprintf_Lg "#include <stdio.h>\nint main(){ long double a=1.0L; char b[128]; snprintf(b,sizeof(b),\"%.40Lg\",a); return 0; }")
endif()

# --- Phase 3: math functions ----------------------------------------------
set(_mp_f128_math 0)
if(_mp_have_f128)
  _mp_try(_mp_f128_math "#include <math.h>\nint main(){ _Float128 x=2.0; (void)sqrtf128(x);(void)expf128(x);(void)logf128(x);(void)sinf128(x);(void)cosf128(x); return 0; }")
endif()

set(_mp_ldbl_math 0)
if(_mp_ldbl_b128)
  _mp_try(_mp_ldbl_math "#include <math.h>\nint main(){ long double x=2.0L; (void)sqrtl(x);(void)expl(x);(void)logl(x);(void)sinl(x);(void)cosl(x); return 0; }")
endif()

# --- Phase 4: literal suffixes --------------------------------------------
set(_mp_suffix_F128 0)
if(_mp_have_f128)
  _mp_try(_mp_suffix_F128 "int main(){ _Float128 a=3.14159265358979323846264338327950288F128; (void)a; return 0; }")
endif()
set(_mp_suffix_Q 0)
if(_mp_have_q128)
  _mp_try(_mp_suffix_Q "int main(){ __float128 a=3.14159265358979323846264338327950288Q; (void)a; return 0; }")
endif()

# --- Phase 5: pick the best configuration (same order as configure.ac) -----
set(_mode 0)
set(_io 0)
set(_math 0)
set(_configured 0)

if(_mp_have_f128 AND _mp_strfromf128 AND _mp_suffix_F128)
  set(_mode 2)   # FLOAT128
  if(_mp_f128_math)
    set(_math 2) # F128
  elseif(_mp_f128_is_ldbl AND _mp_ldbl_math)
    set(_math 1) # LDBL fallback
  else()
    message(FATAL_ERROR "binary128: _Float128 available but no math functions")
  endif()
  set(_io 2)     # STRFROMF128
  set(_configured 1)
  message(STATUS "  binary128: _Float128 + strfromf128")
endif()

if(NOT _configured AND _mp_have_q128 AND _mp_quadmath_snprintf AND _mp_quadmath_math AND _mp_suffix_Q)
  set(_mode 3)   # QUADMATH
  set(_io 3)
  set(_math 3)
  set(_configured 1)
  message(STATUS "  binary128: __float128 + libquadmath")
endif()

if(NOT _configured AND _mp_ldbl_b128 AND _mp_snprintf_Lg AND _mp_ldbl_math)
  set(_mode 1)   # LDBL
  set(_io 1)
  set(_math 1)
  set(_configured 1)
  message(STATUS "  binary128: long double (binary128)")
endif()

if(NOT _configured)
  message(FATAL_ERROR "binary128 is not supported by this compiler/system. "
                      "Use GCC, or pass -DMPLAPACK_ENABLE_BINARY128=OFF.")
endif()

# std::abs availability for the selected binary128 scalar type.
set(_mp_std_abs 0)
if(_mode EQUAL 2 AND _mp_have_f128)
  _mp_try(_mp_std_abs "#include <cmath>
int main(){ auto p=static_cast<_Float128(*)(_Float128)>(&std::abs); (void)p; return 0; }")
elseif(_mode EQUAL 3 AND _mp_have_q128)
  _mp_try(_mp_std_abs "#include <cmath>
int main(){ auto p=static_cast<__float128(*)(__float128)>(&std::abs); (void)p; return 0; }")
endif()

# std::_Float128 scalar math overloads.  These must be probed separately from
# the C f128 functions because modern libstdc++ exposes them by C++ mode.
set(_mp_std_math 0)
if(_mode EQUAL 2)
  _mp_try(_mp_std_math "#include <cmath>
#include <type_traits>
int main(){ _Float128 x=(_Float128)1.0; _Float128 y=(_Float128)2.0; static_assert(std::is_same<decltype(std::sin(x)), _Float128>::value, \"std::sin(_Float128)\"); static_assert(std::is_same<decltype(std::sinh(x)), _Float128>::value, \"std::sinh(_Float128)\"); static_assert(std::is_same<decltype(std::cos(x)), _Float128>::value, \"std::cos(_Float128)\"); static_assert(std::is_same<decltype(std::cosh(x)), _Float128>::value, \"std::cosh(_Float128)\"); static_assert(std::is_same<decltype(std::atan2(x,y)), _Float128>::value, \"std::atan2(_Float128)\"); static_assert(std::is_same<decltype(std::exp(x)), _Float128>::value, \"std::exp(_Float128)\"); static_assert(std::is_same<decltype(std::log(x)), _Float128>::value, \"std::log(_Float128)\"); static_assert(std::is_same<decltype(std::log10(x)), _Float128>::value, \"std::log10(_Float128)\"); static_assert(std::is_same<decltype(std::log2(x)), _Float128>::value, \"std::log2(_Float128)\"); static_assert(std::is_same<decltype(std::pow(x,y)), _Float128>::value, \"std::pow(_Float128)\"); static_assert(std::is_same<decltype(std::sqrt(x)), _Float128>::value, \"std::sqrt(_Float128)\"); static_assert(std::is_same<decltype(std::ceil(x)), _Float128>::value, \"std::ceil(_Float128)\"); static_assert(std::is_same<decltype(std::nextafter(x,y)), _Float128>::value, \"std::nextafter(_Float128)\"); static_assert(std::is_same<decltype(std::ldexp(x,1)), _Float128>::value, \"std::ldexp(_Float128)\"); return 0; }")
endif()

# Complex _Float128 math options used by mplapack_utils_binary128.h.
set(_mp_std_complex 0)
set(_mp_c_complex 0)
if(_mode EQUAL 2)
  _mp_try(_mp_std_complex "#include <complex>
int main(){ std::complex<_Float128> z((_Float128)1.0, (_Float128)2.0); volatile _Float128 a=std::abs(z); (void)a; (void)std::sqrt(z); (void)std::sin(z); (void)std::cos(z); (void)std::exp(z); (void)std::log(z); return 0; }")
  _mp_try(_mp_c_complex "#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif
#include <complex.h>
int main(){ _Float128 _Complex z; __real__(z)=(_Float128)1.0; __imag__(z)=(_Float128)2.0; volatile _Float128 a=cabsf128(z); (void)a; (void)csqrtf128(z); (void)csinf128(z); (void)ccosf128(z); (void)cexpf128(z); (void)clogf128(z); return 0; }")
endif()

# Export results. (This module is include()'d, so it shares the caller's scope;
# plain set() — not PARENT_SCOPE — publishes to the top-level CMakeLists.)
set(MPLAPACK_HAVE__FLOAT128 ${_mp_have_f128})
set(MPLAPACK_HAVE___FLOAT128 ${_mp_have_q128})
set(MPLAPACK_LDBL_IS_BINARY128 ${_mp_ldbl_b128})
set(MPLAPACK__FLOAT128_IS_LDBL ${_mp_f128_is_ldbl})
set(MPLAPACK_FLOAT128_SAME_TYPES ${_mp_same})
set(MPLAPACK_FLOAT128_CONVERT_GNU_TO_STD ${_mp_gnu2std})
set(MPLAPACK_FLOAT128_CONVERT_STD_TO_GNU ${_mp_std2gnu})
set(MPLAPACK_BINARY128_MODE ${_mode})
set(MPLAPACK_BINARY128_IO ${_io})
set(MPLAPACK_BINARY128_MATH ${_math})
set(MPLAPACK_HAVE_STD_ABS_FLOAT128 ${_mp_std_abs})
set(MPLAPACK_HAVE_STD_MATH_FLOAT128 ${_mp_std_math})
set(MPLAPACK_HAVE_STD_COMPLEX_FLOAT128 ${_mp_std_complex})
set(MPLAPACK_HAVE_C_COMPLEX_FLOAT128 ${_mp_c_complex})

# Backend link/compile needs: QUADMATH mode requires -lquadmath.
if(_mode EQUAL 3)
  set(MPLAPACK_BINARY128_EXTRA_LIBS "quadmath")
endif()
