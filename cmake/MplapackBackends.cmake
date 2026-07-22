# MplapackBackends.cmake
#
# Defines the base and optimized self-contained backend libraries. Source lists
# (MPBLAS_SOURCES, MPLAPACK_SOURCES) and shared include directories are set by
# the top-level CMakeLists before this file is included.

function(mplapack_configure_backend_target target macro)
  target_compile_features(${target} PUBLIC cxx_std_17)
  set_target_properties(${target} PROPERTIES
    CXX_STANDARD ${MPLAPACK_CXX_STANDARD}
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS ${MPLAPACK_CXX_EXTENSIONS}
    VERSION ${MPLAPACK_LIB_VERSION}
    SOVERSION ${MPLAPACK_LIB_SOVERSION}
    WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_compile_definitions(${target} PUBLIC ${macro})
  target_include_directories(${target} PUBLIC
    "$<BUILD_INTERFACE:${MPLAPACK_INCLUDE_BUILD_DIRS}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mplapack>")

  if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
  endif()
endfunction()

function(mplapack_add_backend backend macro)
  set(_target mplapack_${backend})
  add_library(${_target} ${MPBLAS_SOURCES} ${MPLAPACK_SOURCES})
  add_library(mplapack::${_target} ALIAS ${_target})
  mplapack_configure_backend_target(${_target} ${macro})
endfunction()

# Remove base sources whose basename is provided by a public override. Helper
# implementations keep suffixes such as _ref and _omp and never shadow.
function(mplapack_shadow_sources out base_list override_list)
  set(_shadow_stems "")
  foreach(_override IN LISTS ${override_list})
    get_filename_component(_stem "${_override}" NAME_WE)
    if(NOT _stem MATCHES "(_ref|_omp|_cuda|_opencl)$")
      list(APPEND _shadow_stems "${_stem}")
    endif()
  endforeach()

  set(_result "")
  foreach(_base IN LISTS ${base_list})
    get_filename_component(_stem "${_base}" NAME_WE)
    if(NOT _stem IN_LIST _shadow_stems)
      list(APPEND _result "${_base}")
    endif()
  endforeach()
  set(${out} "${_result}" PARENT_SCOPE)
endfunction()

function(mplapack_add_opt_backend backend)
  string(TOUPPER "${backend}" _backend_upper)
  set(_macro ___MPLAPACK_BUILD_WITH_${_backend_upper}___)
  set(_target mplapack_${backend}_opt)

  # CPU-only sources. Accelerator subdirectories are separate libraries.
  file(GLOB _mpblas_opt_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${backend}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${backend}/openmp/*.cpp")
  # Despite its helper-style name, this legacy file defines public Rgemm
  # again. The explicit autotools source lists intentionally exclude it.
  list(FILTER _mpblas_opt_sources EXCLUDE REGEX "/openmp/Rgemm_omp\\.cpp$")
  file(GLOB _mplapack_opt_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/mplapack/optimized/${backend}/*.cpp")

  mplapack_shadow_sources(
    _mpblas_reference_sources MPBLAS_SOURCES _mpblas_opt_sources)
  mplapack_shadow_sources(
    _mplapack_reference_sources MPLAPACK_SOURCES _mplapack_opt_sources)

  add_library(${_target}
    ${_mpblas_reference_sources}
    ${_mpblas_opt_sources}
    ${_mplapack_reference_sources}
    ${_mplapack_opt_sources})
  add_library(mplapack::${_target} ALIAS ${_target})
  mplapack_configure_backend_target(${_target} ${_macro})
endfunction()
