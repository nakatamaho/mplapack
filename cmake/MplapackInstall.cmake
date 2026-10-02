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

# mpfrc++ headers are needed by the MPFR backend headers.
if(MPLAPACK_ENABLE_MPFR)
  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/mpfrc++/mpreal.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpfrc++/mpcomplex.h"
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
  string(REGEX REPLACE "^mplapack_" "" _component "${_target}")
  list(APPEND MPLAPACK_AVAILABLE_COMPONENTS "${_component}")
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
  set(PC_LIBS_PRIVATE "")
  set(PC_CFLAGS_EXTRA "")
  set(PC_LIBS_EXTRA "")
  if(PC_FLAVOR MATCHES "^(gmp|mpfr)($|_)")
    foreach(_include IN LISTS GMP_INCLUDE_DIRS GMPXX_INCLUDE_DIR)
      string(APPEND PC_CFLAGS_EXTRA " -I${_include}")
    endforeach()
    get_filename_component(_dep_libdir "${GMP_LIBRARY}" DIRECTORY)
    string(APPEND PC_LIBS_EXTRA " -L${_dep_libdir} -lgmpxx -lgmp")
    if(GMP_PKGCONFIG_FOUND)
      string(JOIN " " PC_LIBS_PRIVATE ${PC_GMP_STATIC_LDFLAGS})
      string(JOIN " " _dependency_options ${PC_GMP_CFLAGS_OTHER})
      string(APPEND PC_CFLAGS_EXTRA " ${_dependency_options}")
    endif()
  endif()
  if(PC_FLAVOR MATCHES "^mpfr($|_)")
    foreach(_include IN LISTS MPC_INCLUDE_DIRS MPFR_INCLUDE_DIRS)
      string(APPEND PC_CFLAGS_EXTRA " -I${_include}")
    endforeach()
    get_filename_component(_mpc_libdir "${MPC_LIBRARY}" DIRECTORY)
    get_filename_component(_mpfr_libdir "${MPFR_LIBRARY}" DIRECTORY)
    set(PC_LIBS_EXTRA "-L${_mpc_libdir} -lmpc -L${_mpfr_libdir} -lmpfr ${PC_LIBS_EXTRA}")
  endif()
  if(PC_FLAVOR MATCHES "^(qd|dd)($|_)")
    foreach(_include IN LISTS QD_INCLUDE_DIRS)
      string(APPEND PC_CFLAGS_EXTRA " -I${_include}")
    endforeach()
    get_filename_component(_dep_libdir "${QD_LIBRARY}" DIRECTORY)
    string(APPEND PC_LIBS_EXTRA " -L${_dep_libdir} -lqd")
    set(PC_LIBS_PRIVATE "-lm")
    if(QD_PKGCONFIG_FOUND)
      string(JOIN " " PC_LIBS_PRIVATE ${PC_QD_STATIC_LDFLAGS})
      string(APPEND PC_LIBS_PRIVATE " -lm")
      string(JOIN " " _dependency_options ${PC_QD_CFLAGS_OTHER})
      string(APPEND PC_CFLAGS_EXTRA " ${_dependency_options}")
    endif()
  endif()
  if(_target MATCHES "_opt($|_)" AND TARGET OpenMP::OpenMP_CXX)
    foreach(_runtime IN LISTS OpenMP_CXX_LIB_NAMES)
      get_filename_component(_runtime_dir "${OpenMP_${_runtime}_LIBRARY}" DIRECTORY)
      if(_runtime_dir)
        string(APPEND PC_LIBS_PRIVATE " -L${_runtime_dir}")
      endif()
      string(APPEND PC_LIBS_PRIVATE " -l${_runtime}")
    endforeach()
    string(STRIP "${PC_LIBS_PRIVATE}" PC_LIBS_PRIVATE)
  endif()
  if(PC_FLAVOR MATCHES "^dd($|_)" AND MPLAPACK_HAS_FP_CONTRACT_OFF)
    string(APPEND PC_CFLAGS_EXTRA " -ffp-contract=off")
  endif()
  string(STRIP "${PC_CFLAGS_EXTRA}" PC_CFLAGS_EXTRA)
  string(STRIP "${PC_LIBS_EXTRA}" PC_LIBS_EXTRA)
  # Do not repeat compiler-default search directories in installed metadata.
  # This matches Autotools' default-system dependency representation.
  foreach(_include IN LISTS CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES)
    string(REPLACE "-I${_include} " "" PC_CFLAGS_EXTRA "${PC_CFLAGS_EXTRA} ")
  endforeach()
  foreach(_directory IN LISTS CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES)
    string(REPLACE "-L${_directory} " "" PC_LIBS_EXTRA "${PC_LIBS_EXTRA} ")
    string(REPLACE "-L${_directory} " "" PC_LIBS_PRIVATE "${PC_LIBS_PRIVATE} ")
  endforeach()
  string(REGEX REPLACE " +" " " PC_CFLAGS_EXTRA "${PC_CFLAGS_EXTRA}")
  string(REGEX REPLACE " +" " " PC_LIBS_EXTRA "${PC_LIBS_EXTRA}")
  string(REGEX REPLACE " +" " " PC_LIBS_PRIVATE "${PC_LIBS_PRIVATE}")
  string(STRIP "${PC_CFLAGS_EXTRA}" PC_CFLAGS_EXTRA)
  string(STRIP "${PC_LIBS_EXTRA}" PC_LIBS_EXTRA)
  string(STRIP "${PC_LIBS_PRIVATE}" PC_LIBS_PRIVATE)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mplapack.pc.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/${_target}.pc"
    @ONLY)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/${_target}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endforeach()
