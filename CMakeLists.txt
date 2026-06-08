cmake_minimum_required(VERSION 3.20)

# -----------------------------------------------------------------------------
# MPLAPACK — multiple-precision BLAS and LAPACK in C++.
#
# This is a modern CMake build that produces, for each enabled precision
# backend, a pair of libraries (mpblas_<backend>, mplapack_<backend>) with
# namespaced targets and an installed package config, so that downstream
# projects can consume MPLAPACK either via find_package(mplapack) or via
# FetchContent.  Precision dependencies (GMP/MPFR/MPC/QD) are located with
# find_package; they are NOT downloaded or built here.
# -----------------------------------------------------------------------------

# Keep version in sync with configure.ac (m4 MPLAPACK_VER_* and LAPACK_VER_*).
project(mplapack
  VERSION 2.2.1
  DESCRIPTION "A multiple precision version of BLAS and LAPACK"
  HOMEPAGE_URL "https://github.com/nakatamaho/mplapack"
  LANGUAGES CXX C)

set(MPLAPACK_LAPACK_VERSION "3.12.1")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Build everything (libraries, tests, examples, benchmarks) at gnu++17, matching
# the autotools build's single -std=gnu++17.  This is the default standard for
# targets in this project's scope; it does NOT leak to a parent project that
# pulls mplapack in via add_subdirectory/FetchContent (such a consumer instead
# inherits the C++17 *floor* exported as a PUBLIC compile feature on the
# backend targets).  Extensions stay on for the __float80/__float128 backends.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# ABI (libtool -version-info 1:0:0) — independent of the package version.
set(MPLAPACK_LIB_VERSION   "1.0.0")
set(MPLAPACK_LIB_SOVERSION "1")

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
option(MPLAPACK_ENABLE_DOUBLE    "Build the double (binary64) backend"      ON)
option(MPLAPACK_ENABLE_MPFR      "Build the MPFR backend"                   ON)
option(MPLAPACK_ENABLE_GMP       "Build the GMP backend"                    ON)
option(MPLAPACK_ENABLE_QD        "Build the QD (quad-double) backend"       ON)
option(MPLAPACK_ENABLE_DD        "Build the DD (double-double) backend"     ON)
option(MPLAPACK_ENABLE_BINARY128 "Build the binary128 backend"             ON)
option(MPLAPACK_ENABLE_BINARY80  "Build the binary80 backend (x86 only)"  OFF)

option(MPLAPACK_BUILD_EXAMPLES   "Build the examples"                      OFF)
option(MPLAPACK_BUILD_TESTS      "Build the test suite"                    OFF)
option(MPLAPACK_BUILD_BENCHMARKS "Build the benchmarks (needs OpenBLAS)"   OFF)

# -----------------------------------------------------------------------------
# Dependencies (find only; never fetch/build)
# -----------------------------------------------------------------------------
# GMP (C++ iface) is needed directly by the GMP backend and transitively by MPFR.
if(MPLAPACK_ENABLE_GMP OR MPLAPACK_ENABLE_MPFR)
  find_package(GMP REQUIRED)
  if(NOT GMPXX_FOUND)
    message(FATAL_ERROR "libgmpxx (the GMP C++ interface) is required for the "
                        "GMP/MPFR backends but was not found.")
  endif()
endif()
if(MPLAPACK_ENABLE_MPFR)
  find_package(MPFR REQUIRED)
  find_package(MPC REQUIRED)
endif()
if(MPLAPACK_ENABLE_QD OR MPLAPACK_ENABLE_DD)
  find_package(QD REQUIRED)
endif()

find_package(OpenMP QUIET COMPONENTS CXX)

# -----------------------------------------------------------------------------
# binary128 / binary80 capability detection -> config-header substitutions
# -----------------------------------------------------------------------------
include(MplapackDetectBinary128)
include(MplapackDetectBinary80)

# -----------------------------------------------------------------------------
# Generated headers: mplapack_config.h and iMlaver.cpp
# -----------------------------------------------------------------------------
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mplapack_config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/mplapack_config.h"
  @ONLY)

set(mplapack_ver_major ${PROJECT_VERSION_MAJOR})
set(mplapack_ver_minor ${PROJECT_VERSION_MINOR})
set(mplapack_ver_patch ${PROJECT_VERSION_PATCH})
string(REPLACE "." ";" _lapack_ver_list "${MPLAPACK_LAPACK_VERSION}")
list(GET _lapack_ver_list 0 LAPACK_VER_MAJOR)
list(GET _lapack_ver_list 1 LAPACK_VER_MINOR)
list(GET _lapack_ver_list 2 LAPACK_VER_PATCH)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mplapack/reference/iMlaver.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/iMlaver.cpp"
  @ONLY)

# Include directories shared by every backend target (BUILD interface).
# The build dir comes first so the CMake-generated mplapack_config.h always
# wins over any stale in-source copy left behind by an autotools configure
# (which writes include/mplapack_config.h into the source tree).
set(MPLAPACK_INCLUDE_BUILD_DIRS
  "${CMAKE_CURRENT_BINARY_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/include")

# -----------------------------------------------------------------------------
# Source lists (shared by all backends)
# -----------------------------------------------------------------------------
file(GLOB MPBLAS_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/reference/*.cpp")

file(GLOB MPLAPACK_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/mplapack/reference/*.cpp")
# iMlaver.cpp is generated from iMlaver.cpp.in; drop any in-tree copy and use ours.
list(FILTER MPLAPACK_SOURCES EXCLUDE REGEX "/iMlaver\\.cpp$")
list(APPEND MPLAPACK_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/generated/iMlaver.cpp")

include(MplapackBackends)

# Track the targets we create so install/export can pick them up.
set(MPLAPACK_INSTALL_TARGETS "")

# A flag every backend can use; only DD strictly needs it, but it is harmless.
check_cxx_compiler_flag("-ffp-contract=off" MPLAPACK_HAS_FP_CONTRACT_OFF)

# -----------------------------------------------------------------------------
# Backends
# -----------------------------------------------------------------------------
if(MPLAPACK_ENABLE_GMP)
  mplapack_add_backend(gmp ___MPLAPACK_BUILD_WITH_GMP___)
  target_link_libraries(mpblas_gmp   PUBLIC GMP::GMPXX GMP::GMP)
  target_link_libraries(mplapack_gmp PUBLIC GMP::GMPXX GMP::GMP)
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_gmp mplapack_gmp)
endif()

if(MPLAPACK_ENABLE_MPFR)
  mplapack_add_backend(mpfr ___MPLAPACK_BUILD_WITH_MPFR___)
  # mpfrc++ headers (mpreal.h / mpcomplex.h) are bundled in-tree.
  foreach(tgt mpblas_mpfr mplapack_mpfr)
    target_include_directories(${tgt} PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpfrc++>")
    target_link_libraries(${tgt} PUBLIC MPC::MPC MPFR::MPFR GMP::GMPXX GMP::GMP)
  endforeach()
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_mpfr mplapack_mpfr)
endif()

if(MPLAPACK_ENABLE_QD)
  mplapack_add_backend(qd ___MPLAPACK_BUILD_WITH_QD___)
  target_link_libraries(mpblas_qd   PUBLIC QD::QD)
  target_link_libraries(mplapack_qd PUBLIC QD::QD)
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_qd mplapack_qd)
endif()

if(MPLAPACK_ENABLE_DD)
  mplapack_add_backend(dd ___MPLAPACK_BUILD_WITH_DD___)
  target_link_libraries(mpblas_dd   PUBLIC QD::QD)
  target_link_libraries(mplapack_dd PUBLIC QD::QD)
  # DD requires FMA contraction to be disabled for correctness.
  if(MPLAPACK_HAS_FP_CONTRACT_OFF)
    target_compile_options(mpblas_dd   PRIVATE -ffp-contract=off)
    target_compile_options(mplapack_dd PRIVATE -ffp-contract=off)
  endif()
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_dd mplapack_dd)
endif()

if(MPLAPACK_ENABLE_DOUBLE)
  mplapack_add_backend(double ___MPLAPACK_BUILD_WITH_DOUBLE___)
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_double mplapack_double)
endif()

if(MPLAPACK_ENABLE_BINARY80)
  mplapack_add_backend(binary80 ___MPLAPACK_BUILD_WITH_BINARY80___)
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_binary80 mplapack_binary80)
endif()

if(MPLAPACK_ENABLE_BINARY128)
  mplapack_add_backend(binary128 ___MPLAPACK_BUILD_WITH_BINARY128___)
  if(MPLAPACK_BINARY128_EXTRA_LIBS)
    target_link_libraries(mpblas_binary128   PUBLIC ${MPLAPACK_BINARY128_EXTRA_LIBS})
    target_link_libraries(mplapack_binary128 PUBLIC ${MPLAPACK_BINARY128_EXTRA_LIBS})
  endif()
  list(APPEND MPLAPACK_INSTALL_TARGETS mpblas_binary128 mplapack_binary128)
endif()

if(NOT MPLAPACK_INSTALL_TARGETS)
  message(FATAL_ERROR "No MPLAPACK backend enabled; enable at least one "
                      "MPLAPACK_ENABLE_<backend> option.")
endif()

# -----------------------------------------------------------------------------
# Optional components
# -----------------------------------------------------------------------------
if(MPLAPACK_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if(MPLAPACK_BUILD_TESTS)
  enable_testing()
  # Shared by both test directories: the override libraries/objects rely on the
  # linker keeping the first (overriding) definition of duplicated symbols.
  include(CheckLinkerFlag)
  check_linker_flag(CXX "-Wl,--allow-multiple-definition" MPLAPACK_LD_ALLOW_MULTIDEF)
  # The mpblas correctness tests compare against a reference Fortran BLAS.
  find_package(BLAS)
  add_subdirectory(mpblas/test)
  add_subdirectory(mplapack/test)
endif()
if(MPLAPACK_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()

# -----------------------------------------------------------------------------
# Install + export
# -----------------------------------------------------------------------------
include(cmake/MplapackInstall.cmake)
