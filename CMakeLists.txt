cmake_minimum_required(VERSION 3.20)

# -----------------------------------------------------------------------------
# MPLAPACK — multiple-precision BLAS and LAPACK in C++.
#
# This is a modern CMake build that produces, for each enabled precision
# backend, one self-contained library (mplapack_<backend>) with
# namespaced targets and an installed package config, so that downstream
# projects can consume MPLAPACK either via find_package(mplapack) or via
# FetchContent.  Precision dependencies (GMP/MPFR/MPC/QD) are located with
# find_package; they are NOT downloaded or built here.
# -----------------------------------------------------------------------------

# Keep version in sync with configure.ac (m4 MPLAPACK_VER_* and LAPACK_VER_*).
project(mplapack
  VERSION 2.3.0
  DESCRIPTION "A multiple precision version of BLAS and LAPACK"
  HOMEPAGE_URL "https://github.com/nakatamaho/mplapack"
  LANGUAGES CXX C)

set(PROJECT_VERSION "${PROJECT_VERSION}-dev")

set(MPLAPACK_LAPACK_VERSION "3.12.1")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Build everything (libraries, tests, examples, benchmarks) with a fixed C++
# standard.  The default remains gnu++17.  Callers can select the standard with
# -DMPLAPACK_CXX_STANDARD=20/23 or -DCMAKE_CXX_STANDARD=20/23, and can select
# strict ISO mode with -DMPLAPACK_CXX_EXTENSIONS=OFF or -DCMAKE_CXX_EXTENSIONS=OFF.
# Consumers still inherit only the C++17 floor exported as a PUBLIC compile
# feature on the backend targets.
if(DEFINED CMAKE_CXX_STANDARD AND NOT CMAKE_CXX_STANDARD STREQUAL "")
  set(_mplapack_default_cxx_standard "${CMAKE_CXX_STANDARD}")
else()
  set(_mplapack_default_cxx_standard 17)
endif()
if(DEFINED CMAKE_CXX_EXTENSIONS AND NOT CMAKE_CXX_EXTENSIONS STREQUAL "")
  set(_mplapack_default_cxx_extensions "${CMAKE_CXX_EXTENSIONS}")
else()
  set(_mplapack_default_cxx_extensions ON)
endif()
set(MPLAPACK_CXX_STANDARD "${_mplapack_default_cxx_standard}" CACHE STRING
    "C++ standard used to build MPLAPACK (17, 20, or 23)")
set_property(CACHE MPLAPACK_CXX_STANDARD PROPERTY STRINGS 17 20 23)
option(MPLAPACK_CXX_EXTENSIONS
       "Build MPLAPACK with compiler C++ extensions enabled (gnu++ mode)"
       ${_mplapack_default_cxx_extensions})
if(NOT MPLAPACK_CXX_STANDARD MATCHES "^(17|20|23)$")
  message(FATAL_ERROR "MPLAPACK_CXX_STANDARD must be 17, 20, or 23")
endif()
set(CMAKE_CXX_STANDARD ${MPLAPACK_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ${MPLAPACK_CXX_EXTENSIONS})

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
option(MPLAPACK_ENABLE_OPT       "Build optimized libraries"                ON)
option(MPLAPACK_ENABLE_CUDA      "Build the DD CUDA accelerator library"    OFF)
option(MPLAPACK_ENABLE_OPENCL    "Build the binary128 OpenCL library"       OFF)
set(MPLAPACK_CUDA_ARCHITECTURES "70;80;90" CACHE STRING
    "CUDA architectures for mplapack_dd_opt_cuda")

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
if(MPLAPACK_ENABLE_CUDA)
  if(NOT MPLAPACK_ENABLE_DD)
    message(FATAL_ERROR "MPLAPACK_ENABLE_CUDA requires MPLAPACK_ENABLE_DD=ON")
  endif()
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # glibc 2.42 exposes C23 rsqrt declarations that conflict with CUDA 13.1.
    # X/Open keeps the POSIX declarations required by the C++ standard library.
    string(APPEND CMAKE_CUDA_FLAGS_INIT
      " -U_GNU_SOURCE -D_XOPEN_SOURCE=700")
  endif()
  enable_language(CUDA)
endif()
if(MPLAPACK_ENABLE_OPENCL)
  if(NOT MPLAPACK_ENABLE_BINARY128)
    message(FATAL_ERROR
      "MPLAPACK_ENABLE_OPENCL requires MPLAPACK_ENABLE_BINARY128=ON")
  endif()
  find_package(OpenCL REQUIRED)
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
  target_link_libraries(mplapack_gmp PUBLIC GMP::GMPXX GMP::GMP)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_gmp)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(gmp)
    target_link_libraries(mplapack_gmp_opt PUBLIC GMP::GMPXX GMP::GMP)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_gmp_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_MPFR)
  mplapack_add_backend(mpfr ___MPLAPACK_BUILD_WITH_MPFR___)
  # mpfrc++ headers (mpreal.h / mpcomplex.h) are bundled in-tree.
  target_include_directories(mplapack_mpfr PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpfrc++>")
  target_link_libraries(mplapack_mpfr PUBLIC MPC::MPC MPFR::MPFR GMP::GMPXX GMP::GMP)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_mpfr)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(mpfr)
    target_include_directories(mplapack_mpfr_opt PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpfrc++>")
    target_link_libraries(mplapack_mpfr_opt PUBLIC
        MPC::MPC MPFR::MPFR GMP::GMPXX GMP::GMP)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_mpfr_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_QD)
  mplapack_add_backend(qd ___MPLAPACK_BUILD_WITH_QD___)
  target_link_libraries(mplapack_qd PUBLIC QD::QD)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_qd)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(qd)
    target_link_libraries(mplapack_qd_opt PUBLIC QD::QD)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_qd_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_DD)
  mplapack_add_backend(dd ___MPLAPACK_BUILD_WITH_DD___)
  target_link_libraries(mplapack_dd PUBLIC QD::QD)
  # DD requires FMA contraction to be disabled for correctness.
  if(MPLAPACK_HAS_FP_CONTRACT_OFF)
    target_compile_options(mplapack_dd PRIVATE -ffp-contract=off)
  endif()
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_dd)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(dd)
    target_link_libraries(mplapack_dd_opt PUBLIC QD::QD)
    if(MPLAPACK_HAS_FP_CONTRACT_OFF)
      target_compile_options(mplapack_dd_opt PRIVATE -ffp-contract=off)
    endif()
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_dd_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_CUDA)
  file(GLOB _mplapack_cuda_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/*.cu")
  # Split kernels are textually included by the six CUDA translation units.
  list(FILTER _mplapack_cuda_sources EXCLUDE REGEX
    "/(Rgemm_(fermi|tesla)_(NN|NT|TN|TT)_[0p]|Rsyrk_[NT][LU]_[0p])\\.cu$")
  # cuda/Rgemm.cu and cuda/Rsyrk.cu shadow the same-stem CPU/reference files.
  mplapack_add_accelerator_backend(dd cuda ${_mplapack_cuda_sources})
  target_sources(mplapack_dd_opt_cuda PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/dd_real_cuda.h")
  target_include_directories(mplapack_dd_opt_cuda PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda")
  target_link_libraries(mplapack_dd_opt_cuda PUBLIC QD::QD CUDA::cudart)
  if(MPLAPACK_HAS_FP_CONTRACT_OFF)
    target_compile_options(mplapack_dd_opt_cuda PRIVATE
      $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
  endif()
  set_target_properties(mplapack_dd_opt_cuda PROPERTIES
    CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON)
  mplapack_add_unique_symbol_post_build(mplapack_dd_opt_cuda)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_dd_opt_cuda)
endif()

if(MPLAPACK_ENABLE_DOUBLE)
  mplapack_add_backend(double ___MPLAPACK_BUILD_WITH_DOUBLE___)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_double)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(double)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_double_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_BINARY80)
  mplapack_add_backend(binary80 ___MPLAPACK_BUILD_WITH_BINARY80___)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary80)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(binary80)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary80_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_BINARY128)
  mplapack_add_backend(binary128 ___MPLAPACK_BUILD_WITH_BINARY128___)
  if(MPLAPACK_BINARY128_EXTRA_LIBS)
    target_link_libraries(mplapack_binary128 PUBLIC ${MPLAPACK_BINARY128_EXTRA_LIBS})
  endif()
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary128)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(binary128)
    if(MPLAPACK_BINARY128_EXTRA_LIBS)
      target_link_libraries(mplapack_binary128_opt PUBLIC
        ${MPLAPACK_BINARY128_EXTRA_LIBS})
    endif()
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary128_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_OPENCL)
  file(GLOB _mplapack_opencl_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/binary128/opencl/*.cpp")
  set_source_files_properties(${_mplapack_opencl_sources} PROPERTIES
    COMPILE_DEFINITIONS CL_TARGET_OPENCL_VERSION=120)
  mplapack_add_accelerator_backend(
    binary128 opencl ${_mplapack_opencl_sources})
  target_sources(mplapack_binary128_opt_opencl PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/binary128/opencl/OclHelper.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/binary128/opencl/cpp_kernel.h")
  target_include_directories(mplapack_binary128_opt_opencl PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/binary128/opencl")
  target_link_libraries(mplapack_binary128_opt_opencl PUBLIC OpenCL::OpenCL)
  if(MPLAPACK_BINARY128_EXTRA_LIBS)
    target_link_libraries(mplapack_binary128_opt_opencl PUBLIC
      ${MPLAPACK_BINARY128_EXTRA_LIBS})
  endif()
  mplapack_add_unique_symbol_post_build(mplapack_binary128_opt_opencl)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary128_opt_opencl)
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
  foreach(b gmp mpfr qd dd double binary80 binary128)
    if(TARGET mplapack_${b})
      add_test(NAME unique_symbols_${b}
        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
                "$<TARGET_FILE:mplapack_${b}>")
    endif()
    if(TARGET mplapack_${b}_opt)
      add_test(NAME unique_symbols_${b}_opt
        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
                "$<TARGET_FILE:mplapack_${b}_opt>")
    endif()
  endforeach()
  foreach(accelerator dd_opt_cuda binary128_opt_opencl)
    if(TARGET mplapack_${accelerator})
      add_test(NAME unique_symbols_${accelerator}
        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
                "$<TARGET_FILE:mplapack_${accelerator}>")
    endif()
  endforeach()
  # Shared by both test directories: the override libraries/objects rely on the
  # linker keeping the first (overriding) definition of duplicated symbols.
  include(CheckLinkerFlag)
  check_linker_flag(CXX "-Wl,--allow-multiple-definition" MPLAPACK_LD_ALLOW_MULTIDEF)
  # The mpblas/mplapack correctness tests compare against a reference Fortran
  # BLAS/LAPACK and need its Fortran runtime.
  find_package(BLAS)
  find_package(LAPACK)
  find_library(MPLAPACK_GFORTRAN_LIBRARY NAMES gfortran)
  add_subdirectory(mpblas/test)
  add_subdirectory(mplapack/test)
  add_subdirectory(mplapack/test/compare)
endif()
if(MPLAPACK_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()

# -----------------------------------------------------------------------------
# Install + export
# -----------------------------------------------------------------------------
include(cmake/MplapackInstall.cmake)

# -----------------------------------------------------------------------------
# Packaging (source tarball; the `make dist` equivalent) — only when this is the
# top-level project, so a parent FetchContent build does not inherit our CPack.
# -----------------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(cmake/MplapackPackage.cmake)
endif()
