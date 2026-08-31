# MplapackInstall.cmake — install rules, target export, package config, pkg-config.
# Included from the top-level CMakeLists after all backend targets are defined.

include(CMakePackageConfigHelpers)

set(MPLAPACK_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/mplapack"
    CACHE STRING "Install location for mplapack CMake package files")

# --- Libraries -------------------------------------------------------------
install(TARGETS ${MPLAPACK_INSTALL_TARGETS}
  EXPORT mplapackTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR})

# --- Headers ---------------------------------------------------------------
# Match the autotools install layout: public headers live in
# <prefix>/include/mplapack and are included as, for example,
# <mplapack_mpfr.h> through the exported target's include directory.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mplapack"
  FILES_MATCHING PATTERN "*.h")

# Generated config header.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/mplapack_config.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mplapack")

if(MPLAPACK_ENABLE_GMP OR MPLAPACK_ENABLE_MPFR)
  install(DIRECTORY "${MPLAPACK_GMPFRXX_MKII_ROOT}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mplapack")
endif()

# --- Export targets --------------------------------------------------------
install(EXPORT mplapackTargets
  FILE mplapackTargets.cmake
  NAMESPACE mplapack::
  DESTINATION "${MPLAPACK_INSTALL_CMAKEDIR}")

# Make the build tree usable via FetchContent without installing.
export(EXPORT mplapackTargets
  NAMESPACE mplapack::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/mplapackTargets.cmake")

# --- Package config --------------------------------------------------------
# Which dependencies the consumer must re-find.
if(MPLAPACK_ENABLE_GMP OR MPLAPACK_ENABLE_MPFR)
  set(MPLAPACK_NEEDS_GMP 1)
else()
  set(MPLAPACK_NEEDS_GMP 0)
endif()
if(MPLAPACK_ENABLE_MPFR)
  set(MPLAPACK_NEEDS_MPFR 1)
else()
  set(MPLAPACK_NEEDS_MPFR 0)
endif()
if(MPLAPACK_ENABLE_QD OR MPLAPACK_ENABLE_DD)
  set(MPLAPACK_NEEDS_QD 1)
else()
  set(MPLAPACK_NEEDS_QD 0)
endif()
if(TARGET OpenMP::OpenMP_CXX)
  set(MPLAPACK_NEEDS_OPENMP 1)
else()
  set(MPLAPACK_NEEDS_OPENMP 0)
endif()
if(TARGET mplapack_dd_opt_cuda)
  set(MPLAPACK_NEEDS_CUDA 1)
else()
  set(MPLAPACK_NEEDS_CUDA 0)
endif()
if(TARGET mplapack_binary128_opt_opencl)
  set(MPLAPACK_NEEDS_OPENCL 1)
else()
  set(MPLAPACK_NEEDS_OPENCL 0)
endif()

set(MPLAPACK_ENABLED_BACKENDS "")
foreach(b gmp mpfr qd dd double binary80 binary128)
  string(TOUPPER ${b} B)
  if(MPLAPACK_ENABLE_${B})
    list(APPEND MPLAPACK_ENABLED_BACKENDS ${b})
  endif()
endforeach()

set(MPLAPACK_AVAILABLE_COMPONENTS "")
foreach(_target IN LISTS MPLAPACK_INSTALL_TARGETS)
  if(_target MATCHES "^mplapack_")
    string(REGEX REPLACE "^mplapack_" "" _component "${_target}")
    list(APPEND MPLAPACK_AVAILABLE_COMPONENTS "${_component}")
  endif()
endforeach()

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mplapackConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mplapackConfig.cmake"
  INSTALL_DESTINATION "${MPLAPACK_INSTALL_CMAKEDIR}")

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/mplapackConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/mplapackConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/mplapackConfigVersion.cmake"
  DESTINATION "${MPLAPACK_INSTALL_CMAKEDIR}")

# Ship the bundled Find modules so find_dependency() works on the consumer side.
install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGMP.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindMPFR.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindMPC.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindQD.cmake"
  DESTINATION "${MPLAPACK_INSTALL_CMAKEDIR}")

# The build-tree package config uses the same dependency discovery path as the
# installed package, so copy the bundled Find modules next to it as well.
file(COPY
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGMP.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindMPFR.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindMPC.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindQD.cmake"
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# --- pkg-config (per-flavor only; identical convention in autotools) -------
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pkgconfig")
foreach(_target IN LISTS MPLAPACK_INSTALL_TARGETS)
  set(PC_PREFIX "${CMAKE_INSTALL_PREFIX}")
  set(PC_LIBDIR "${CMAKE_INSTALL_FULL_LIBDIR}")
  set(PC_INCLUDEDIR "${CMAKE_INSTALL_FULL_INCLUDEDIR}/mplapack")
  set(PC_NAME "${_target}")
  string(REGEX REPLACE "^mplapack_" "" PC_FLAVOR "${_target}")
  set(PC_DESCRIPTION "${PROJECT_DESCRIPTION}")
  set(PC_VERSION "${PROJECT_VERSION}")
  if(_target STREQUAL "mplapack_mpfr" OR
     _target STREQUAL "mplapack_mpfr_opt")
    set(PC_REQUIRES "mpc mpfr gmp")
  else()
    set(PC_REQUIRES "")
  endif()
  if(_target STREQUAL "mplapack_mpfr_opt" AND OpenMP_CXX_FOUND)
    set(PC_LIBS_PRIVATE "${OpenMP_CXX_FLAGS}")
  else()
    set(PC_LIBS_PRIVATE "")
  endif()
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mplapack.pc.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/${_target}.pc"
    @ONLY)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/${_target}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endforeach()
