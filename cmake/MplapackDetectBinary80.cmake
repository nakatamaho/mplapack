# MplapackDetectBinary80.cmake
#
# Reproduces the binary80 (80-bit x87 extended) probing from configure.ac.
# Sets, in the caller's scope:
#
#   MPLAPACK_HAVE__FLOAT64X            MPLAPACK_HAVE___FLOAT80
#   MPLAPACK_LDBL_IS_BINARY80          MPLAPACK_FLOAT64X_IS_BINARY80
#   MPLAPACK__FLOAT64X_IS_LDBL         MPLAPACK_FLOAT80_SAME_TYPES
#   MPLAPACK_FLOAT80_CONVERT_GNU_TO_STD  MPLAPACK_FLOAT80_CONVERT_STD_TO_GNU
#   MPLAPACK_BINARY80_MODE  (0 DISABLED / 1 LDBL80 / 2 FLOAT64X)
#   MPLAPACK_BINARY80_IO    (0 NONE / 1 SNPRINTF_LDBL / 2 STRFROMF64X)
#   MPLAPACK_BINARY80_MATH  (0 NONE / 1 LDBL / 2 F64X)
#
# Priority follows configure.ac: _Float64x(binary80) > long double(binary80).

include(CheckCXXSourceCompiles)

function(_mp80_try var src)
  check_cxx_source_compiles("${src}" ${var}_RESULT)
  if(${var}_RESULT)
    set(${var} 1 PARENT_SCOPE)
  else()
    set(${var} 0 PARENT_SCOPE)
  endif()
endfunction()

set(MPLAPACK_HAVE__FLOAT64X 0)
set(MPLAPACK_HAVE___FLOAT80 0)
set(MPLAPACK_LDBL_IS_BINARY80 0)
set(MPLAPACK_FLOAT64X_IS_BINARY80 0)
set(MPLAPACK__FLOAT64X_IS_LDBL 0)
set(MPLAPACK_FLOAT80_SAME_TYPES 0)
set(MPLAPACK_FLOAT80_CONVERT_GNU_TO_STD 0)
set(MPLAPACK_FLOAT80_CONVERT_STD_TO_GNU 0)
set(MPLAPACK_BINARY80_MODE 0)
set(MPLAPACK_BINARY80_IO 0)
set(MPLAPACK_BINARY80_MATH 0)

if(NOT MPLAPACK_ENABLE_BINARY80)
  return()
endif()

message(STATUS "Detecting binary80 support")

# --- Phase 1: types --------------------------------------------------------
_mp80_try(_b80_f64x "int main(){ _Float64x a,b,c; b=1.0; c=3.0; a=b/c; (void)a; return 0; }")
_mp80_try(_b80_f80  "int main(){ __float80 a,b,c; b=1.0; c=3.0; a=b/c; (void)a; return 0; }")
_mp80_try(_b80_ldbl "int main(){\n#if (__LDBL_MANT_DIG__ == 64) && (__LDBL_MAX_EXP__ == 16384)\nreturn 0;\n#else\n#error not binary80\n#endif\n}")

set(_b80_f64x_is_b80 0)
if(_b80_f64x)
  _mp80_try(_b80_f64x_is_b80 "int main(){\n#if defined(__FLT64X_MANT_DIG__) && defined(__FLT64X_MAX_EXP__)\n#if (__FLT64X_MANT_DIG__ == 64) && (__FLT64X_MAX_EXP__ == 16384)\nreturn 0;\n#else\n#error\n#endif\n#else\n#error\n#endif\n}")
endif()

set(_b80_f64x_is_ldbl 0)
if(_b80_f64x)
  _mp80_try(_b80_f64x_is_ldbl "#include <type_traits>\nint main(){ static_assert(std::is_same<_Float64x,long double>::value,\"\"); return 0; }")
endif()

set(_b80_same 0)
set(_b80_gnu2std 0)
set(_b80_std2gnu 0)
if(_b80_f64x AND _b80_f80)
  _mp80_try(_b80_same "#include <type_traits>\nint main(){ static_assert(std::is_same<_Float64x,__float80>::value,\"\"); return 0; }")
  _mp80_try(_b80_gnu2std "#include <type_traits>\nint main(){ static_assert(std::is_convertible<__float80,_Float64x>::value,\"\"); return 0; }")
  _mp80_try(_b80_std2gnu "#include <type_traits>\nint main(){ static_assert(std::is_convertible<_Float64x,__float80>::value,\"\"); return 0; }")
endif()

# --- Phase 2: I/O ----------------------------------------------------------
set(_b80_strfromf64x 0)
if(_b80_f64x AND _b80_f64x_is_b80)
  _mp80_try(_b80_strfromf64x "#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1\n#define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1\n#include <stdlib.h>\nint main(){ _Float64x a=1.0; char b[128]; strfromf64x(b,sizeof(b),\"%.40g\",a); return 0; }")
endif()

set(_b80_snprintf_Lg 0)
if(_b80_ldbl)
  _mp80_try(_b80_snprintf_Lg "#include <stdio.h>\nint main(){ long double a=1.0L; char b[128]; snprintf(b,sizeof(b),\"%.40Lg\",a); return 0; }")
endif()

# --- Phase 3: math ---------------------------------------------------------
set(_b80_f64x_math 0)
if(_b80_f64x AND _b80_f64x_is_b80)
  _mp80_try(_b80_f64x_math "#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1\n#define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1\n#include <math.h>\nint main(){ _Float64x x=2.0; (void)sqrtf64x(x);(void)expf64x(x);(void)logf64x(x);(void)sinf64x(x);(void)cosf64x(x); return 0; }")
endif()

set(_b80_ldbl_math 0)
if(_b80_ldbl)
  _mp80_try(_b80_ldbl_math "#include <math.h>\nint main(){ long double x=2.0L; (void)sqrtl(x);(void)expl(x);(void)logl(x);(void)sinl(x);(void)cosl(x); return 0; }")
endif()

# --- Phase 5: choose -------------------------------------------------------
set(_mode 0)
set(_io 0)
set(_math 0)
set(_configured 0)

if(_b80_f64x AND _b80_f64x_is_b80 AND _b80_strfromf64x AND _b80_f64x_math)
  set(_mode 2)
  set(_io 2)
  set(_math 2)
  set(_configured 1)
  message(STATUS "  binary80: _Float64x (binary80)")
endif()

if(NOT _configured AND _b80_ldbl AND _b80_snprintf_Lg AND _b80_ldbl_math)
  set(_mode 1)
  set(_io 1)
  set(_math 1)
  set(_configured 1)
  message(STATUS "  binary80: long double (binary80)")
endif()

if(NOT _configured)
  message(FATAL_ERROR "binary80 is not supported by this compiler/system (x86 only). "
                      "Pass -DMPLAPACK_ENABLE_BINARY80=OFF.")
endif()

# Export (include()'d module shares caller scope; plain set()).
set(MPLAPACK_HAVE__FLOAT64X ${_b80_f64x})
set(MPLAPACK_HAVE___FLOAT80 ${_b80_f80})
set(MPLAPACK_LDBL_IS_BINARY80 ${_b80_ldbl})
set(MPLAPACK_FLOAT64X_IS_BINARY80 ${_b80_f64x_is_b80})
set(MPLAPACK__FLOAT64X_IS_LDBL ${_b80_f64x_is_ldbl})
set(MPLAPACK_FLOAT80_SAME_TYPES ${_b80_same})
set(MPLAPACK_FLOAT80_CONVERT_GNU_TO_STD ${_b80_gnu2std})
set(MPLAPACK_FLOAT80_CONVERT_STD_TO_GNU ${_b80_std2gnu})
set(MPLAPACK_BINARY80_MODE ${_mode})
set(MPLAPACK_BINARY80_IO ${_io})
set(MPLAPACK_BINARY80_MATH ${_math})
