# Benchmarks.  Each R<op> driver is compiled once per enabled backend and linked
# against the matching mplapack_<backend> library (which carries the
# ___MPLAPACK_BUILD_WITH_<backend>___ definition and the MPBLAS layer).  The
# drivers dlopen() an optional optimized library at run time, hence ${CMAKE_DL_LIBS}.
#
# The double-precision reference/OpenBLAS comparison programs from the autotools
# build are intentionally out of scope here (they need a Fortran reference BLAS
# and a bundled OpenBLAS).

set(_bench_drivers Raxpy Rcopy Rdot Rgemm Rgemv Rsyrk Rpotrf Rpotf2 Rgetrf)
set(_bench_backends gmp mpfr qd dd double binary80 binary128)

foreach(_b ${_bench_backends})
  string(TOUPPER ${_b} _B)
  if(NOT MPLAPACK_ENABLE_${_B})
    continue()
  endif()
  foreach(_drv ${_bench_drivers})
    set(_src "${CMAKE_CURRENT_SOURCE_DIR}/${_drv}.cpp")
    if(NOT EXISTS "${_src}")
      continue()
    endif()
    set(_name "bench_${_drv}_${_b}")
    add_executable(${_name} "${_src}")
    target_link_libraries(${_name} PRIVATE mplapack::mplapack_${_b} ${CMAKE_DL_LIBS})
    if(TARGET OpenMP::OpenMP_CXX)
      target_link_libraries(${_name} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set_target_properties(${_name} PROPERTIES
      OUTPUT_NAME "${_drv}.${_b}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endforeach()
endforeach()
