# MPBLAS correctness tests, ported from the autotools mpblas/test harness.
#
# Each common/<Routine>.test.cpp compares one MPBLAS routine against the
# reference Fortran BLAS. Test support always uses MPFR for random-number
# generation, so these tests require the MPFR backend and a reference BLAS.

if(NOT MPLAPACK_ENABLE_MPFR)
  message(STATUS "mpblas/test: skipped (needs the MPFR backend for its test support)")
  return()
endif()
if(NOT BLAS_FOUND)
  message(STATUS "mpblas/test: skipped (no reference BLAS found)")
  return()
endif()

set(_common "${CMAKE_CURRENT_SOURCE_DIR}/common")
set(_support
  "${_common}/xerbla.cpp"
  "${_common}/mplapack.test.cpp"
  "${_common}/Mxerbla.override.cpp")

file(GLOB _mpblas_tests CONFIGURE_DEPENDS "${_common}/*.test.cpp")

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

    foreach(_src ${_mpblas_tests})
      get_filename_component(_stem "${_src}" NAME_WE)
      set(_exe "${_stem}_test_${_variant}")

      add_executable(${_exe} "${_src}" ${_support})
      # Set one backend macro explicitly. Linking target files avoids inheriting
      # the MPFR oracle target's conflicting backend macro.
      target_compile_definitions(${_exe} PRIVATE
        ___MPLAPACK_BUILD_WITH_${_B}___ ___MPLAPACK_INTERNAL___)
      target_include_directories(${_exe} PRIVATE
        "${CMAKE_BINARY_DIR}/include"
        "${CMAKE_SOURCE_DIR}/include"
        "${CMAKE_SOURCE_DIR}/mpfrc++")

      set(_backend_archives "$<TARGET_FILE:${_backend_library}>")
      set(_backend_dependencies ${_backend_library})
      if(NOT _b STREQUAL "mpfr")
        list(APPEND _backend_archives "$<TARGET_FILE:mplapack_mpfr>")
        list(APPEND _backend_dependencies mplapack_mpfr)
      endif()
      target_link_libraries(${_exe} PRIVATE
        ${_backend_archives}
        MPC::MPC MPFR::MPFR GMP::GMPXX
        ${BLAS_LIBRARIES})
      if(_b STREQUAL "dd" OR _b STREQUAL "qd")
        target_link_libraries(${_exe} PRIVATE QD::QD)
      endif()
      if(_b STREQUAL "binary128" AND MPLAPACK_BINARY128_EXTRA_LIBS)
        target_link_libraries(${_exe} PRIVATE ${MPLAPACK_BINARY128_EXTRA_LIBS})
      endif()
      add_dependencies(${_exe} ${_backend_dependencies})

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

      add_test(NAME mpblas_${_variant}_${_stem} COMMAND ${_exe})
    endforeach()
  endforeach()
endforeach()
