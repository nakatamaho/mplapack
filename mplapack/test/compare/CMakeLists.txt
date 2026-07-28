# MPLAPACK "compare" tests, ported from the autotools mplapack/test/compare
# harness.  Each common/<Routine>.test.cpp exercises one MPLAPACK routine and
# compares it against the reference Fortran LAPACK/BLAS (dgesv_ etc., declared in
# lapack.h/blas.h via F77_FUNC); it exit(1)s on a mismatch.
#
# Like mpblas/test, the support code always uses MPFR, so this directory needs
# the MPFR backend and reference BLAS + LAPACK; it is skipped otherwise.

if(NOT MPLAPACK_ENABLE_MPFR)
  message(STATUS "mplapack/test/compare: skipped (needs the MPFR backend)")
  return()
endif()
if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(STATUS "mplapack/test/compare: skipped (no reference BLAS/LAPACK found)")
  return()
endif()

set(_common  "${CMAKE_CURRENT_SOURCE_DIR}/common")
set(_mptest_common "${CMAKE_SOURCE_DIR}/mpblas/test/common")

# Shared support (same files as mpblas/test): xerbla_ for the reference library,
# the MPFR-based random/compare helpers, and the override Mxerbla which wins from
# the executable's own objects via --allow-multiple-definition.
set(_support
  "${_mptest_common}/xerbla.cpp"
  "${_mptest_common}/mplapack.test.cpp"
  "${_mptest_common}/Mxerbla.override.cpp")

file(GLOB _compare_tests CONFIGURE_DEPENDS "${_common}/*.test.cpp")

set(_backends gmp mpfr qd dd double binary80 binary128)
foreach(_b ${_backends})
  string(TOUPPER ${_b} _B)
  if(NOT MPLAPACK_ENABLE_${_B})
    continue()
  endif()

  set(_backend_libraries mplapack_${_b})
  if(TARGET mplapack_${_b}_opt)
    list(APPEND _backend_libraries mplapack_${_b}_opt)
  endif()

  foreach(_backend_library IN LISTS _backend_libraries)
    if(_backend_library MATCHES "_opt$")
      set(_variant "${_b}_opt")
    else()
      set(_variant "${_b}")
    endif()

    foreach(_src ${_compare_tests})
      get_filename_component(_stem "${_src}" NAME_WE)
      set(_exe "${_stem}_compare_${_variant}")

      set(_extra "")
      if(_stem STREQUAL "Mutils")
        set(_extra "${_common}/ilaenv.cpp")   # Mutils.test needs the test ILAENV
      endif()

      add_executable(${_exe} "${_src}" ${_support} ${_extra})
      # Exactly one backend macro (the merged archives are linked by file path
      # below so their target macros do not leak in).
      target_compile_definitions(${_exe} PRIVATE
        MPLAPACK_BUILD_WITH_${_B} MPLAPACK_INTERNAL)
      if(_b STREQUAL "dd" AND MPLAPACK_HAS_FP_CONTRACT_OFF)
        target_compile_options(${_exe} PRIVATE -ffp-contract=off)
      endif()
      target_include_directories(${_exe} PRIVATE
        "${CMAKE_BINARY_DIR}/include"
        "${CMAKE_SOURCE_DIR}/include"
        "${MPLAPACK_GMPFRXX_MKII_ROOT}/include")

      set(_backend_archives "$<TARGET_FILE:${_backend_library}>")
      set(_backend_dependencies ${_backend_library})
      if(NOT _b STREQUAL "mpfr")
        list(APPEND _backend_archives "$<TARGET_FILE:mplapack_mpfr>")
        list(APPEND _backend_dependencies mplapack_mpfr)
      endif()
      target_link_libraries(${_exe} PRIVATE
        ${_backend_archives}
        MPC::MPC MPFR::MPFR         ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
      add_dependencies(${_exe} ${_backend_dependencies})
      if(_b STREQUAL "dd" OR _b STREQUAL "qd")
        target_link_libraries(${_exe} PRIVATE QD::QD)
      endif()
      if(_b STREQUAL "binary128" AND MPLAPACK_BINARY128_EXTRA_LIBS)
        target_link_libraries(${_exe} PRIVATE ${MPLAPACK_BINARY128_EXTRA_LIBS})
      endif()
      if(MPLAPACK_GFORTRAN_LIBRARY)
        target_link_libraries(${_exe} PRIVATE ${MPLAPACK_GFORTRAN_LIBRARY})
      endif()
      if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(${_exe} PRIVATE OpenMP::OpenMP_CXX)
      endif()
      if(MPLAPACK_LD_ALLOW_MULTIDEF)
        target_link_options(${_exe} PRIVATE "-Wl,--allow-multiple-definition")
      endif()
      set_target_properties(${_exe} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_variant}")

      add_test(NAME compare_${_variant}_${_stem} COMMAND ${_exe})
    endforeach()
  endforeach()
endforeach()
