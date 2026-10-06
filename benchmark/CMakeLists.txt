# Benchmarks.  Each R<op> driver is compiled once per enabled backend and linked
# against the matching mplapack_<backend> library (which carries the
# MPLAPACK_BUILD_WITH_<backend> definition and the MPBLAS layer), and again
# against mplapack_<backend>_opt and the CUDA libraries when they are built.  The
# drivers dlopen() an optional optimized library at run time, hence ${CMAKE_DL_LIBS}.

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

# Optimized libraries: <op>.<backend>_opt (as the autotools build).
foreach(_b ${_bench_backends})
  if(NOT TARGET mplapack_${_b}_opt)
    continue()
  endif()
  foreach(_drv ${_bench_drivers})
    set(_src "${CMAKE_CURRENT_SOURCE_DIR}/${_drv}.cpp")
    if(NOT EXISTS "${_src}")
      continue()
    endif()
    set(_name "bench_${_drv}_${_b}_opt")
    add_executable(${_name} "${_src}")
    target_link_libraries(${_name} PRIVATE mplapack::mplapack_${_b}_opt ${CMAKE_DL_LIBS})
    if(TARGET OpenMP::OpenMP_CXX)
      target_link_libraries(${_name} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set_target_properties(${_name} PROPERTIES
      OUTPUT_NAME "${_drv}.${_b}_opt"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endforeach()
endforeach()

# CUDA libraries: <op>.<backend>_cuda_total (host-to-host timings, Rgemm on
# the GPU; algorithms chosen at run time, see go.Rgemm_cuda.sh.in).
foreach(_b dd qd mpfr)
  if(NOT TARGET mplapack_${_b}_opt_cuda)
    continue()
  endif()
  foreach(_drv Rgemm Rpotrf Rgetrf)
    set(_name "bench_${_drv}_${_b}_cuda_total")
    add_executable(${_name} "${CMAKE_CURRENT_SOURCE_DIR}/${_drv}.cpp")
    target_link_libraries(${_name} PRIVATE mplapack::mplapack_${_b}_opt_cuda ${CMAKE_DL_LIBS})
    if(TARGET OpenMP::OpenMP_CXX)
      target_link_libraries(${_name} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set_target_properties(${_name} PROPERTIES
      OUTPUT_NAME "${_drv}.${_b}_cuda_total"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endforeach()
endforeach()

# Runner scripts go.<op>.sh from go.<op>.sh.in, with the gnuplot templates
# next to them (the autotools build does the same at "make" time).
string(TOLOWER "${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}" _bench_cc_tag)
file(GLOB _bench_go_templates "${CMAKE_CURRENT_SOURCE_DIR}/go.*.sh.in")
foreach(_in ${_bench_go_templates})
  get_filename_component(_fn "${_in}" NAME)
  string(REGEX REPLACE "\\.in$" "" _out "${_fn}")
  file(READ "${_in}" _text)
  string(REPLACE "%%CC_TAG%%" "${_bench_cc_tag}" _text "${_text}")
  string(REPLACE "%%CPU_OS_DETECTION_SH%%" "${CMAKE_SOURCE_DIR}/misc/cpu_os_detection.sh" _text "${_text}")
  string(REPLACE "%%PREFIX%%" "${CMAKE_INSTALL_PREFIX}" _text "${_text}")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${_out}" "${_text}")
  file(CHMOD "${CMAKE_CURRENT_BINARY_DIR}/${_out}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_in}")
endforeach()
file(GLOB _bench_plt_templates "${CMAKE_CURRENT_SOURCE_DIR}/*.plt.in")
file(COPY ${_bench_plt_templates} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# -----------------------------------------------------------------------------
# Double-precision comparison benchmarks: each d<op>.cpp is timed against the
# reference Fortran BLAS/LAPACK (find_package) as d<op>_ref, and against OpenBLAS
# as d<op>_openblas.  These do not use any mplapack backend.  Each variant is
# built only when its library is available.
# -----------------------------------------------------------------------------
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_library(MPLAPACK_OPENBLAS_LIBRARY NAMES openblas)
find_library(MPLAPACK_BENCH_GFORTRAN NAMES gfortran)

set(_blas_bench   daxpy dcopy ddot dgemm dgemv dsyrk)  # need reference BLAS
set(_lapack_bench dgetrf dpotf2 dpotrf)                # need reference LAPACK

function(_mplapack_add_dbench op needs_lapack)
  set(_src "${CMAKE_CURRENT_SOURCE_DIR}/${op}.cpp")
  if(NOT EXISTS "${_src}")
    return()
  endif()

  # Reference BLAS/LAPACK variant.
  set(_ref_ok ${BLAS_FOUND})
  if(needs_lapack AND NOT LAPACK_FOUND)
    set(_ref_ok FALSE)
  endif()
  if(_ref_ok)
    add_executable(${op}_ref "${_src}")
    target_include_directories(${op}_ref PRIVATE
      "${CMAKE_BINARY_DIR}/include" "${CMAKE_SOURCE_DIR}/include")
    target_link_libraries(${op}_ref PRIVATE ${CMAKE_DL_LIBS})
    if(needs_lapack)
      target_link_libraries(${op}_ref PRIVATE ${LAPACK_LIBRARIES})
    endif()
    target_link_libraries(${op}_ref PRIVATE ${BLAS_LIBRARIES})
    set_target_properties(${op}_ref PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  # OpenBLAS variant.
  if(MPLAPACK_OPENBLAS_LIBRARY)
    add_executable(${op}_openblas "${_src}")
    target_include_directories(${op}_openblas PRIVATE
      "${CMAKE_BINARY_DIR}/include" "${CMAKE_SOURCE_DIR}/include")
    target_link_libraries(${op}_openblas PRIVATE
      ${MPLAPACK_OPENBLAS_LIBRARY} ${CMAKE_DL_LIBS})
    set_target_properties(${op}_openblas PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  foreach(_v ref openblas)
    if(TARGET ${op}_${_v})
      if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(${op}_${_v} PRIVATE OpenMP::OpenMP_CXX)
      endif()
      if(MPLAPACK_BENCH_GFORTRAN)
        target_link_libraries(${op}_${_v} PRIVATE ${MPLAPACK_BENCH_GFORTRAN})
      endif()
    endif()
  endforeach()
endfunction()

foreach(_op ${_blas_bench})
  _mplapack_add_dbench(${_op} FALSE)
endforeach()
foreach(_op ${_lapack_bench})
  _mplapack_add_dbench(${_op} TRUE)
endforeach()
