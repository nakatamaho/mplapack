# MplapackBackends.cmake
#
# Defines mplapack_add_backend(), which creates one self-contained library for
# a precision backend:
#
#   mplapack_<backend>   (MPBLAS + MPLAPACK sources)
#                        -> alias mplapack::mplapack_<backend>
#
# The target carries the PUBLIC compile definition that selects the backend
# (___MPLAPACK_BUILD_WITH_<BACKEND>___) so that consumers including <mplapack.h>
# automatically resolve to the matching backend header.
#
# Source lists (MPBLAS_SOURCES, MPLAPACK_SOURCES) and the include directories
# are expected to be set by the top-level CMakeLists before this is used.

function(mplapack_add_backend backend macro)
  set(mplapack_tgt mplapack_${backend})

  add_library(${mplapack_tgt} ${MPBLAS_SOURCES} ${MPLAPACK_SOURCES})
  add_library(mplapack::${mplapack_tgt} ALIAS ${mplapack_tgt})

  # PUBLIC floor: consumers compiling against the mplapack headers need at
  # least C++17 (the headers use std::enable_if_t / is_same_v / fold exprs).
  target_compile_features(${mplapack_tgt} PUBLIC cxx_std_17)
  # Pin this library compilation to the project-selected standard. The
  # cxx_std_17 floor alone permits newer compiler defaults; this property keeps
  # library builds reproducible while MPLAPACK_CXX_EXTENSIONS controls whether
  # targets use gnu++XX or strict c++XX mode.
  set_target_properties(${mplapack_tgt} PROPERTIES
    CXX_STANDARD ${MPLAPACK_CXX_STANDARD}
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS ${MPLAPACK_CXX_EXTENSIONS})
  # The backend-selecting macro must reach consumers -> PUBLIC.
  target_compile_definitions(${mplapack_tgt} PUBLIC ${macro})
  target_include_directories(${mplapack_tgt} PUBLIC
    "$<BUILD_INTERFACE:${MPLAPACK_INCLUDE_BUILD_DIRS}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mplapack>")
  set_target_properties(${mplapack_tgt} PROPERTIES
    VERSION ${MPLAPACK_LIB_VERSION}
    SOVERSION ${MPLAPACK_LIB_SOVERSION}
    WINDOWS_EXPORT_ALL_SYMBOLS ON)

  # MPLAPACK uses OpenMP in a few routines; link it when available.
  if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(${mplapack_tgt} PRIVATE OpenMP::OpenMP_CXX)
  endif()
endfunction()
