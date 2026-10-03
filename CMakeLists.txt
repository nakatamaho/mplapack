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
  VERSION 3.0.1
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
set(MPLAPACK_CHECK_UNIQUE_SYMBOLS_ENV
    "${CMAKE_COMMAND}" "-E" "env"
    "CXX=${CMAKE_CXX_COMPILER}"
    "NM=${CMAKE_NM}"
    "OBJDUMP=${CMAKE_OBJDUMP}")

# ABI (libtool -version-info 3:0:0) — independent of the package version.
set(MPLAPACK_LIB_VERSION   "3.0.1")
set(MPLAPACK_LIB_SOVERSION "3")

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
option(MPLAPACK_ENABLE_MPFR_CUDA "Build the MPFR CUDA accelerator library (Rgemm on the GPU)" OFF)
set(MPLAPACK_MPC_CUDA_INCLUDE_DIR "" CACHE PATH
    "Include directory of an installed mpc_cuda (contains mpc_cuda.cuh); empty: build external/mpc_cuda")
set(MPLAPACK_MPC_CUDA_LIBRARY "" CACHE FILEPATH
    "libmpc_cuda.a of an installed mpc_cuda (default: <include dir>/../lib/libmpc_cuda.a)")
option(MPLAPACK_ENABLE_OPENCL    "Build the binary128 OpenCL library"       OFF)
option(MPLAPACK_ENABLE_GMPFRXX_MKII "Use bundled gmpfrxx_mkII for GMP/MPFR backends" ON)
set(MPLAPACK_CUDA_ARCHITECTURES "70;80;90" CACHE STRING
    "CUDA architectures for mplapack_dd_opt_cuda")

option(MPLAPACK_BUILD_EXAMPLES   "Build the examples"                      OFF)
option(MPLAPACK_BUILD_TESTS      "Build the test suite"                    OFF)
option(MPLAPACK_BUILD_BENCHMARKS "Build the benchmarks (needs OpenBLAS)"   OFF)

set(MPLAPACK_INTEGER_BITS "auto" CACHE STRING
    "Bit width for mplapackint (auto, 32, or 64)")
set_property(CACHE MPLAPACK_INTEGER_BITS PROPERTY STRINGS auto 32 64)
if(MPLAPACK_INTEGER_BITS STREQUAL "auto")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(MPLAPACK_INTEGER_BITS_CONFIG 64)
  else()
    set(MPLAPACK_INTEGER_BITS_CONFIG 32)
  endif()
elseif(MPLAPACK_INTEGER_BITS STREQUAL "32")
  set(MPLAPACK_INTEGER_BITS_CONFIG 32)
elseif(MPLAPACK_INTEGER_BITS STREQUAL "64")
  set(MPLAPACK_INTEGER_BITS_CONFIG 64)
else()
  message(FATAL_ERROR "MPLAPACK_INTEGER_BITS must be auto, 32, or 64")
endif()
if(MPLAPACK_INTEGER_BITS_CONFIG EQUAL 64)
  set(MPLAPACK_INTEGER_TYPE int64_t)
else()
  set(MPLAPACK_INTEGER_TYPE int32_t)
endif()

# -----------------------------------------------------------------------------
# Dependencies and bundled gmpfrxx_mkII header extraction
# -----------------------------------------------------------------------------
# GMP C interface is needed directly by the GMP backend and transitively by MPFR.
if(MPLAPACK_ENABLE_GMP OR MPLAPACK_ENABLE_MPFR)
  find_package(GMP REQUIRED)
endif()
if(MPLAPACK_ENABLE_MPFR)
  find_package(MPFR REQUIRED)
  find_package(MPC REQUIRED)
endif()
if(MPLAPACK_ENABLE_GMP OR MPLAPACK_ENABLE_MPFR)
  if(NOT MPLAPACK_ENABLE_GMPFRXX_MKII)
    message(FATAL_ERROR "GMP/MPFR backends require MPLAPACK_ENABLE_GMPFRXX_MKII=ON")
  endif()
  set(MPLAPACK_GMPFRXX_MKII_ROOT
      "${CMAKE_CURRENT_BINARY_DIR}/external/i/GMPFRXX_MKII")
  set(MPLAPACK_GMPFRXX_MKII_VERSION "1.4.3")
  set(MPLAPACK_GMPFRXX_MKII_ARCHIVE
      "${CMAKE_CURRENT_SOURCE_DIR}/external/gmpfrxx_mkII/download/gmpfrxx_mkII.${MPLAPACK_GMPFRXX_MKII_VERSION}.tar.xz")
  set(MPLAPACK_GMPFRXX_MKII_STAMP
      "${MPLAPACK_GMPFRXX_MKII_ROOT}/.gmpfrxx_mkII-${MPLAPACK_GMPFRXX_MKII_VERSION}")
  if(NOT EXISTS "${MPLAPACK_GMPFRXX_MKII_STAMP}")
    find_program(MPLAPACK_TAR_EXECUTABLE tar REQUIRED)
    file(REMOVE_RECURSE "${MPLAPACK_GMPFRXX_MKII_ROOT}")
    file(MAKE_DIRECTORY "${MPLAPACK_GMPFRXX_MKII_ROOT}")
    execute_process(
      COMMAND "${MPLAPACK_TAR_EXECUTABLE}" -xJf "${MPLAPACK_GMPFRXX_MKII_ARCHIVE}"
              -C "${MPLAPACK_GMPFRXX_MKII_ROOT}" --strip-components=1
      COMMAND_ERROR_IS_FATAL ANY)
    file(WRITE "${MPLAPACK_GMPFRXX_MKII_STAMP}" "${MPLAPACK_GMPFRXX_MKII_VERSION}\n")
  endif()
endif()
if(MPLAPACK_ENABLE_QD OR MPLAPACK_ENABLE_DD)
  find_package(QD REQUIRED)
endif()
if(MPLAPACK_ENABLE_CUDA AND NOT MPLAPACK_ENABLE_DD)
  message(FATAL_ERROR "MPLAPACK_ENABLE_CUDA requires MPLAPACK_ENABLE_DD=ON")
endif()
if(MPLAPACK_ENABLE_MPFR_CUDA AND NOT (MPLAPACK_ENABLE_MPFR AND MPLAPACK_ENABLE_OPT))
  message(FATAL_ERROR
    "MPLAPACK_ENABLE_MPFR_CUDA requires MPLAPACK_ENABLE_MPFR=ON and MPLAPACK_ENABLE_OPT=ON")
endif()
set(_mplapack_mpc_cuda_bundled OFF)
if(MPLAPACK_ENABLE_MPFR_CUDA)
  if(MPLAPACK_MPC_CUDA_INCLUDE_DIR)
    if(NOT EXISTS "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}/mpc_cuda.cuh")
      message(FATAL_ERROR
        "MPLAPACK_MPC_CUDA_INCLUDE_DIR must point to the mpc_cuda include directory "
        "(containing mpc_cuda.cuh); see https://github.com/tkouya/mpc_cuda")
    endif()
    if(NOT MPLAPACK_MPC_CUDA_LIBRARY)
      set(MPLAPACK_MPC_CUDA_LIBRARY "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}/../lib/libmpc_cuda.a")
    endif()
    if(NOT EXISTS "${MPLAPACK_MPC_CUDA_LIBRARY}")
      message(FATAL_ERROR
        "MPLAPACK_MPC_CUDA_LIBRARY must point to libmpc_cuda.a built for MPLAPACK_CUDA_ARCHITECTURES")
    endif()
  else()
    # Build the bundled mpc_cuda (external/mpc_cuda); it targets one architecture.
    list(LENGTH MPLAPACK_CUDA_ARCHITECTURES _mplapack_cuda_arch_count)
    if(NOT _mplapack_cuda_arch_count EQUAL 1)
      message(FATAL_ERROR
        "The bundled mpc_cuda is built for one GPU architecture: set "
        "MPLAPACK_CUDA_ARCHITECTURES to a single value (e.g. 90 or 121), or give an "
        "installed mpc_cuda with MPLAPACK_MPC_CUDA_INCLUDE_DIR")
    endif()
    set(_mplapack_mpc_cuda_bundled ON)
    set(_mplapack_mpc_cuda_root "${CMAKE_CURRENT_BINARY_DIR}/external/i/MPC_CUDA")
    set(MPLAPACK_MPC_CUDA_INCLUDE_DIR "${_mplapack_mpc_cuda_root}/include")
    set(MPLAPACK_MPC_CUDA_LIBRARY "${_mplapack_mpc_cuda_root}/lib/libmpc_cuda.a")
    file(MAKE_DIRECTORY "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
  endif()
endif()
if(MPLAPACK_ENABLE_CUDA OR MPLAPACK_ENABLE_MPFR_CUDA)
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
  mplapack_add_backend(gmp MPLAPACK_BUILD_WITH_GMP)
  target_include_directories(mplapack_gmp PUBLIC
      "$<BUILD_INTERFACE:${MPLAPACK_GMPFRXX_MKII_ROOT}/include>")
  target_link_libraries(mplapack_gmp PUBLIC GMP::GMP)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_gmp)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(gmp)
    target_include_directories(mplapack_gmp_opt PUBLIC
        "$<BUILD_INTERFACE:${MPLAPACK_GMPFRXX_MKII_ROOT}/include>")
    target_link_libraries(mplapack_gmp_opt PUBLIC GMP::GMP)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_gmp_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_MPFR)
  mplapack_add_backend(mpfr MPLAPACK_BUILD_WITH_MPFR)
  target_include_directories(mplapack_mpfr PUBLIC
      "$<BUILD_INTERFACE:${MPLAPACK_GMPFRXX_MKII_ROOT}/include>")
  target_link_libraries(mplapack_mpfr PUBLIC MPC::MPC MPFR::MPFR GMP::GMP)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_mpfr)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(mpfr)
    target_include_directories(mplapack_mpfr_opt PUBLIC
        "$<BUILD_INTERFACE:${MPLAPACK_GMPFRXX_MKII_ROOT}/include>")
    target_link_libraries(mplapack_mpfr_opt PUBLIC MPC::MPC MPFR::MPFR GMP::GMP)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_mpfr_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_QD)
  mplapack_add_backend(qd MPLAPACK_BUILD_WITH_QD)
  target_link_libraries(mplapack_qd PUBLIC QD::QD)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_qd)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(qd)
    target_link_libraries(mplapack_qd_opt PUBLIC QD::QD)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_qd_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_DD)
  mplapack_add_backend(dd MPLAPACK_BUILD_WITH_DD)
  target_link_libraries(mplapack_dd PUBLIC QD::QD)
  mplapack_configure_dd_arithmetic(mplapack_dd)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_dd)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(dd)
    target_link_libraries(mplapack_dd_opt PUBLIC QD::QD)
    mplapack_configure_dd_arithmetic(mplapack_dd_opt)
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
  mplapack_configure_dd_arithmetic(mplapack_dd_opt_cuda)
  set_target_properties(mplapack_dd_opt_cuda PROPERTIES
    CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON)
  mplapack_add_unique_symbol_post_build(mplapack_dd_opt_cuda)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_dd_opt_cuda)
endif()

if(_mplapack_mpc_cuda_bundled)
  include(ExternalProject)
  set(_mplapack_external "${CMAKE_CURRENT_SOURCE_DIR}/external")
  set(_mplapack_mpc_cuda_nvcc "${CMAKE_CUDA_COMPILER}")
  if(CMAKE_CUDA_HOST_COMPILER)
    string(APPEND _mplapack_mpc_cuda_nvcc " -ccbin ${CMAKE_CUDA_HOST_COMPILER}")
  endif()
  # mpc_cuda builds from the GMP/MPFR/MPC sources of external/{gmp,mpfr,mpc}
  # and must not be built in parallel.
  ExternalProject_Add(mpc_cuda_bundled
    URL "${_mplapack_external}/mpc_cuda/download/mpc_cuda-0.0.1-abd9407.tar.xz"
    URL_HASH SHA256=bcac28655a62d381ffe38fc09704aca39f45610d1f697c351b4fc0acaa2278e7
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/external/mpc_cuda"
    PATCH_COMMAND
      ${CMAKE_COMMAND} -E tar xf "${_mplapack_external}/gmp/download/gmp-6.3.0.tar.xz"
      COMMAND ${CMAKE_COMMAND} -E tar xf "${_mplapack_external}/mpfr/download/mpfr-4.2.2.tar.bz2"
      COMMAND ${CMAKE_COMMAND} -E tar xf "${_mplapack_external}/mpc/download/mpc-1.4.1.tar.xz"
    CONFIGURE_COMMAND <SOURCE_DIR>/configure "--prefix=${_mplapack_mpc_cuda_root}"
      "--with-cuda-arch=sm_${MPLAPACK_CUDA_ARCHITECTURES}" "NVCC=${_mplapack_mpc_cuda_nvcc}"
    BUILD_IN_SOURCE ON
    BUILD_COMMAND ${CMAKE_COMMAND} -E env MAKEFLAGS= make -j1
    INSTALL_COMMAND ${CMAKE_COMMAND} -E env MAKEFLAGS= make -j1 install
    BUILD_BYPRODUCTS "${MPLAPACK_MPC_CUDA_LIBRARY}")
endif()

if(MPLAPACK_ENABLE_MPFR_CUDA)
  # cuda/Rgemm.cpp shadows optimized/mpfr/Rgemm.cpp; the *_cuda files are helpers.
  set(_mplapack_mpfr_cuda_dir "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/mpfr/cuda")
  mplapack_add_accelerator_backend(mpfr cuda
    "${_mplapack_mpfr_cuda_dir}/Rgemm.cpp"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_bridge_cuda.cpp"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_device_cuda.cu"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu")
  target_sources(mplapack_mpfr_opt_cuda PRIVATE
    "${_mplapack_mpfr_cuda_dir}/Rgemm_kernel_cuda.h"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.h")
  target_include_directories(mplapack_mpfr_opt_cuda
    PUBLIC "$<BUILD_INTERFACE:${MPLAPACK_GMPFRXX_MKII_ROOT}/include>"
    PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}" "${_mplapack_mpfr_cuda_dir}")
  target_link_libraries(mplapack_mpfr_opt_cuda
    PUBLIC MPC::MPC MPFR::MPFR GMP::GMP CUDA::cudart
    PRIVATE "${MPLAPACK_MPC_CUDA_LIBRARY}")
  # cu_mpfr (Rgemm_rt_cuda.cu) is relocatable device code from libmpc_cuda.a;
  # resolve it inside the library so that C++ programs can link it normally.
  target_compile_options(mplapack_mpfr_opt_cuda PRIVATE
    "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
  set_target_properties(mplapack_mpfr_opt_cuda PROPERTIES
    CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  if(_mplapack_mpc_cuda_bundled)
    add_dependencies(mplapack_mpfr_opt_cuda mpc_cuda_bundled)
  endif()
  mplapack_add_unique_symbol_post_build(mplapack_mpfr_opt_cuda)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_mpfr_opt_cuda)
endif()

if(MPLAPACK_ENABLE_DOUBLE)
  mplapack_add_backend(double MPLAPACK_BUILD_WITH_DOUBLE)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_double)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(double)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_double_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_BINARY80)
  mplapack_add_backend(binary80 MPLAPACK_BUILD_WITH_BINARY80)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary80)
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(binary80)
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_binary80_opt)
  endif()
endif()

if(MPLAPACK_ENABLE_BINARY128)
  mplapack_add_backend(binary128 MPLAPACK_BUILD_WITH_BINARY128)
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
        COMMAND ${MPLAPACK_CHECK_UNIQUE_SYMBOLS_ENV}
                "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
                "$<TARGET_FILE:mplapack_${b}>")
    endif()
    if(TARGET mplapack_${b}_opt)
      add_test(NAME unique_symbols_${b}_opt
        COMMAND ${MPLAPACK_CHECK_UNIQUE_SYMBOLS_ENV}
                "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
                "$<TARGET_FILE:mplapack_${b}_opt>")
    endif()
  endforeach()
  foreach(accelerator dd_opt_cuda mpfr_opt_cuda binary128_opt_opencl)
    if(TARGET mplapack_${accelerator})
      add_test(NAME unique_symbols_${accelerator}
        COMMAND ${MPLAPACK_CHECK_UNIQUE_SYMBOLS_ENV}
                "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
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
  if(MPLAPACK_ENABLE_MPFR_CUDA)
    # Rgemm_cuda_test.cpp checks the GPU Rgemm bit for bit against the CPU
    # Rgemm.  The host variant runs the kernels' element code on the host
    # (cu_mpfr mapped to the system MPFR) and needs no GPU; the device variant
    # is skipped when no device is present.
    # Each runs twice: default dispatch, and with every precision forced to
    # the runtime-precision (cu_mpfr) kernels.  OMP_NUM_THREADS=1 because the
    # OpenMP CPU Rgemm used as the reference creates its temporaries at the
    # worker threads' MPFR default precision.
    set(_mplapack_mpfr_cuda_dir "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/mpfr/cuda")
    foreach(_variant host device)
      if(_variant STREQUAL "host")
        set(_fixed "${_mplapack_mpfr_cuda_dir}/Rgemm_host_cuda.cpp")
      else()
        set(_fixed "${_mplapack_mpfr_cuda_dir}/Rgemm_device_cuda.cu")
      endif()
      set(_exe mpfr_cuda_Rgemm_${_variant}_test)
      add_executable(${_exe}
        "${_mplapack_mpfr_cuda_dir}/test/Rgemm_cuda_test.cpp"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_bridge_cuda.cpp"
        "${_fixed}"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu")
      if(_variant STREQUAL "host")
        target_compile_definitions(${_exe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
      endif()
      target_include_directories(${_exe} PRIVATE
        "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}" "${_mplapack_mpfr_cuda_dir}")
      target_compile_options(${_exe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
      if(_variant STREQUAL "host")
        target_link_libraries(${_exe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      else()
        target_link_libraries(${_exe} PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
      endif()
      set_target_properties(${_exe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON)
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_exe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_Rgemm_${_variant} COMMAND ${_exe})
      add_test(NAME mpfr_cuda_Rgemm_${_variant}_runtime COMMAND ${_exe})
      set_tests_properties(mpfr_cuda_Rgemm_${_variant} PROPERTIES
        SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1")
      set_tests_properties(mpfr_cuda_Rgemm_${_variant}_runtime PROPERTIES
        SKIP_RETURN_CODE 77
        ENVIRONMENT "OMP_NUM_THREADS=1;MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1")

      # Winograd (Strassen) path: exact on integers, bounded on reals, and
      # the GPU backend equal to the host one.
      set(_wexe mpfr_cuda_Rgemm_winograd_${_variant}_test)
      add_executable(${_wexe}
        "${_mplapack_mpfr_cuda_dir}/test/Rgemm_winograd_cuda_test.cpp"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_bridge_cuda.cpp"
        "${_fixed}"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu")
      if(_variant STREQUAL "host")
        target_compile_definitions(${_wexe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
        target_link_libraries(${_wexe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      else()
        target_link_libraries(${_wexe} PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
      endif()
      target_include_directories(${_wexe} PRIVATE
        "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}" "${_mplapack_mpfr_cuda_dir}")
      target_compile_options(${_wexe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
      set_target_properties(${_wexe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON)
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_wexe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_Rgemm_winograd_${_variant} COMMAND ${_wexe})
      set_tests_properties(mpfr_cuda_Rgemm_winograd_${_variant} PROPERTIES
        SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1")
    endforeach()
  endif()
  add_subdirectory(mplapack/test)
  add_subdirectory(mplapack/test/compare)
  add_subdirectory(cmake/tests)
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
