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

# Intel compilers default to a value-unsafe floating-point model (fp-model
# fast: reassociation, contraction), which breaks double-double and the
# other expansion types.  As configure.ac does for CFLAGS/CXXFLAGS, build
# every translation unit with the strict model.
if(CMAKE_CXX_COMPILER_ID MATCHES "^Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fp-model=strict>")
endif()
if(CMAKE_C_COMPILER_ID MATCHES "^Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-fp-model=strict>")
endif()
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
# One MPLAPACK_ENABLE_<NAME> option per row of backends.txt.
include(MplapackBackendTable)
foreach(_b IN LISTS MPLAPACK_BACKENDS)
  string(TOUPPER "${_b}" _B)
  option(MPLAPACK_ENABLE_${_B} "Build the ${_b} backend"
    ${MPLAPACK_BACKEND_${_b}_DEFAULT})
endforeach()
option(MPLAPACK_ENABLE_OPT       "Build optimized libraries"                ON)
option(MPLAPACK_ENABLE_CUDA      "Build the DD CUDA accelerator library"    OFF)
option(MPLAPACK_ENABLE_MPFR_CUDA "Build the MPFR CUDA accelerator library (Rgemm on the GPU)" OFF)
option(MPLAPACK_MPFR_CUDA_RT_LARGE "MPFR CUDA: also the 8192/16384-bit runtime-precision complex BLAS (long device compilation)" OFF)
option(MPLAPACK_ENABLE_QD_CUDA   "Build the QD CUDA accelerator library (Rgemm on the GPU)" OFF)
option(MPLAPACK_ENABLE_TD_CUDA   "Build the TD CUDA accelerator library (Rgemm on the GPU)" OFF)
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
# Dependencies follow the traits of the enabled backends (backends.txt).
mplapack_enabled_backend_traits(_mplapack_enabled_traits)
foreach(_trait gmp mpfr qd)
  string(TOUPPER "${_trait}" _TRAIT)
  if(_trait IN_LIST _mplapack_enabled_traits)
    set(MPLAPACK_USES_${_TRAIT} ON)
  else()
    set(MPLAPACK_USES_${_TRAIT} OFF)
  endif()
endforeach()

# GMP C interface is needed directly by the GMP backend and transitively by MPFR.
if(MPLAPACK_USES_GMP OR MPLAPACK_USES_MPFR)
  find_package(GMP REQUIRED)
endif()
if(MPLAPACK_USES_MPFR)
  find_package(MPFR REQUIRED)
  find_package(MPC REQUIRED)
endif()
if(MPLAPACK_USES_GMP OR MPLAPACK_USES_MPFR)
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
# Keep in step with the libQD3 tarball bundled under external/qd.
set(MPLAPACK_QD_MIN_VERSION 1.6.0)
if(MPLAPACK_USES_QD)
  find_package(QD ${MPLAPACK_QD_MIN_VERSION} REQUIRED)
  # dd, td and qd are computed with libQD3's branch-free arithmetic (QD_BF;
  # mpblas/optimized/common/Rgemm_qd_ops_common.h).
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES "${QD_INCLUDE_DIR}")
  check_cxx_source_compiles([=[
#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#if !defined(QD_BF_ADD) || !defined(QD_BF_MUL)
#error libQD3 without QD_BF
#endif
int main() { return 0; }
]=] MPLAPACK_QD_HAS_BF)
  unset(CMAKE_REQUIRED_INCLUDES)
  if(NOT MPLAPACK_QD_HAS_BF)
    message(FATAL_ERROR
      "libQD3 in ${QD_INCLUDE_DIR} is not built with branch-free arithmetic. "
      "MPLAPACK needs libQD3 configured with -DQD_ENABLE_BF=ON (QD_BF).")
  endif()
endif()
if(MPLAPACK_ENABLE_CUDA AND NOT MPLAPACK_ENABLE_DD)
  message(FATAL_ERROR "MPLAPACK_ENABLE_CUDA requires MPLAPACK_ENABLE_DD=ON")
endif()
if(MPLAPACK_ENABLE_QD_CUDA AND NOT MPLAPACK_ENABLE_QD)
  message(FATAL_ERROR "MPLAPACK_ENABLE_QD_CUDA requires MPLAPACK_ENABLE_QD=ON")
endif()
if(MPLAPACK_ENABLE_TD_CUDA AND NOT MPLAPACK_ENABLE_TD)
  message(FATAL_ERROR "MPLAPACK_ENABLE_TD_CUDA requires MPLAPACK_ENABLE_TD=ON")
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
    if(NOT EXISTS "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}/mpc_cuda/cu_fround.cuh")
      message(FATAL_ERROR "Installed mpc_cuda needs cu_fround.cuh; apply external/mpc_cuda/patches/patch-cu_freal-correctly-rounded-div-sqrt")
    endif()
    if(NOT MPLAPACK_MPC_CUDA_LIBRARY)
      set(MPLAPACK_MPC_CUDA_LIBRARY "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}/../lib/libmpc_cuda.a")
    endif()
    if(NOT EXISTS "${MPLAPACK_MPC_CUDA_LIBRARY}")
      message(FATAL_ERROR
        "MPLAPACK_MPC_CUDA_LIBRARY must point to libmpc_cuda.a built for MPLAPACK_CUDA_ARCHITECTURES")
    endif()
  else()
    # Build the bundled mpc_cuda (external/mpc_cuda) for every architecture of
    # MPLAPACK_CUDA_ARCHITECTURES (e.g. "80;90" for A100 and H100): the first
    # is its --with-cuda-arch, the others are added to its NVCC as -gencode.
    foreach(_arch IN LISTS MPLAPACK_CUDA_ARCHITECTURES)
      if(NOT _arch MATCHES "^[0-9]+[a-z]?$")
        message(FATAL_ERROR
          "The bundled mpc_cuda needs plain compute capabilities in "
          "MPLAPACK_CUDA_ARCHITECTURES (e.g. 80;90), not '${_arch}'")
      endif()
    endforeach()
    set(_mplapack_mpc_cuda_bundled ON)
    set(_mplapack_mpc_cuda_root "${CMAKE_CURRENT_BINARY_DIR}/external/i/MPC_CUDA")
    set(MPLAPACK_MPC_CUDA_INCLUDE_DIR "${_mplapack_mpc_cuda_root}/include")
    set(MPLAPACK_MPC_CUDA_LIBRARY "${_mplapack_mpc_cuda_root}/lib/libmpc_cuda.a")
    file(MAKE_DIRECTORY "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
  endif()
endif()
if(MPLAPACK_ENABLE_CUDA OR MPLAPACK_ENABLE_QD_CUDA OR MPLAPACK_ENABLE_TD_CUDA OR MPLAPACK_ENABLE_MPFR_CUDA)
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
# One library (and one optimized library) per enabled row of backends.txt.
foreach(_b IN LISTS MPLAPACK_BACKENDS)
  string(TOUPPER "${_b}" _B)
  if(NOT MPLAPACK_ENABLE_${_B})
    continue()
  endif()
  mplapack_add_backend(${_b} MPLAPACK_BUILD_WITH_${_B})
  mplapack_configure_backend_traits(mplapack_${_b} ${_b})
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_${_b})
  if(MPLAPACK_ENABLE_OPT)
    mplapack_add_opt_backend(${_b})
    mplapack_configure_backend_traits(mplapack_${_b}_opt ${_b})
    list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_${_b}_opt)
  endif()
endforeach()
# The quad-double kernels of the optimized Rgemm (mpblas/optimized/common)
# rely on uncontracted floating-point operations, as dd and td do.
if(TARGET mplapack_qd_opt)
  mplapack_configure_dd_arithmetic(mplapack_qd_opt)
endif()

# dd/td/qd Level 1/2/3 BLAS on the GPU (dispatchers and host side).
file(GLOB _mplapack_blas_cuda_sources CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/common/blas_cuda/*.cpp")

foreach(_b dd td qd)
  set_property(SOURCE "mpblas/optimized/${_b}/cuda/Rpotrf_gpu_cuda.cu"
    "mpblas/optimized/${_b}/cuda/Rgetrf_gpu_cuda.cu"
    "mpblas/optimized/${_b}/cuda/Rsolver_gpu_cuda.cu" APPEND PROPERTY
    COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-ffp-contract=off>")
endforeach()

foreach(_b dd td qd mpfr)
  set(_solver_${_b}_sources "")
  foreach(_r Rgetrs Rpotrs Rgesv Rposv Cgetrs Cpotrs Cgesv Cposv)
    list(APPEND _solver_${_b}_sources "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${_b}/cuda/${_r}.cpp")
  endforeach()
endforeach()

if(MPLAPACK_ENABLE_CUDA)
  file(GLOB _mplapack_cuda_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/*.cu")
  # Split kernels are textually included by the six CUDA translation units.
  list(FILTER _mplapack_cuda_sources EXCLUDE REGEX
    "/(Rgemm_(fermi|tesla)_(NN|NT|TN|TT)_[0p]|Rsyrk_[NT][LU]_[0p])\\.cu$")
  # cuda/Rgemm.cpp and the Level 1/2/3 BLAS of optimized/common/blas_cuda
  # shadow the same-stem CPU files; Rgemm_gpu_cuda.cu (kernels of
  # mpblas/optimized/common), Rsyrk_legacy_cuda.cu and the *_bridge_cuda.cpp
  # files are helpers.
  mplapack_add_accelerator_backend(dd cuda ${_mplapack_cuda_sources} ${_solver_dd_sources}
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/Rgemm.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/Rgemm_bridge_cuda.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/Rpotrf.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/Cpotrf.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/Rgetrf.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/Cgetrf.cpp"
    ${_mplapack_blas_cuda_sources})
  target_sources(mplapack_dd_opt_cuda PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda/dd_real_cuda.h")
  target_include_directories(mplapack_dd_opt_cuda PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/dd/cuda")
  target_link_libraries(mplapack_dd_opt_cuda PUBLIC QD::QD CUDA::cudart CUDA::cublas)
  mplapack_configure_dd_arithmetic(mplapack_dd_opt_cuda)
  # The double-double device arithmetic must not be contracted into FMAs.
  target_compile_options(mplapack_dd_opt_cuda PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
  set_target_properties(mplapack_dd_opt_cuda PROPERTIES
    CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON)
  mplapack_add_unique_symbol_post_build(mplapack_dd_opt_cuda)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_dd_opt_cuda)
endif()

# qd and td CUDA libraries: Rgemm and Level 1/2/3 BLAS of
# mpblas/optimized/common on the GPU.
foreach(_b qd td)
  string(TOUPPER "${_b}" _B)
  if(NOT MPLAPACK_ENABLE_${_B}_CUDA)
    continue()
  endif()
  # cuda/Rgemm.cpp shadows optimized/<b>/Rgemm.cpp and common/blas_cuda/R*.cpp
  # the CPU dispatchers of common/blas; the *_cuda files are helpers.
  set(_dir "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${_b}/cuda")
  mplapack_add_accelerator_backend(${_b} cuda ${_solver_${_b}_sources} "${_dir}/Rsolver_gpu_cuda.cu"
    "${_dir}/Rgemm.cpp" "${_dir}/Rgemm_bridge_cuda.cpp" "${_dir}/Rgemm_gpu_cuda.cu"
    "${_dir}/Rpotrf.cpp" "${_dir}/Cpotrf.cpp" "${_dir}/Rpotrf_gpu_cuda.cu"
    "${_dir}/Rgetrf.cpp" "${_dir}/Cgetrf.cpp" "${_dir}/Rgetrf_gpu_cuda.cu"
    ${_mplapack_blas_cuda_sources})
  target_link_libraries(mplapack_${_b}_opt_cuda PUBLIC QD::QD CUDA::cudart CUDA::cublas)
  mplapack_configure_dd_arithmetic(mplapack_${_b}_opt_cuda)
  target_compile_options(mplapack_${_b}_opt_cuda PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
  set_target_properties(mplapack_${_b}_opt_cuda PROPERTIES
    CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON)
  mplapack_add_unique_symbol_post_build(mplapack_${_b}_opt_cuda)
  list(APPEND MPLAPACK_INSTALL_TARGETS mplapack_${_b}_opt_cuda)
endforeach()

if(_mplapack_mpc_cuda_bundled)
  include(ExternalProject)
  set(_mplapack_external "${CMAKE_CURRENT_SOURCE_DIR}/external")
  set(_mplapack_mpc_cuda_nvcc "${CMAKE_CUDA_COMPILER}")
  if(CMAKE_CUDA_HOST_COMPILER)
    string(APPEND _mplapack_mpc_cuda_nvcc " -ccbin ${CMAKE_CUDA_HOST_COMPILER}")
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # the glibc 2.42 / CUDA 13.1 rsqrt conflict, as for CMAKE_CUDA_FLAGS_INIT
    string(APPEND _mplapack_mpc_cuda_nvcc " -U_GNU_SOURCE -D_XOPEN_SOURCE=700")
  endif()
  list(GET MPLAPACK_CUDA_ARCHITECTURES 0 _mplapack_mpc_cuda_arch)
  set(_mplapack_mpc_cuda_extra_archs ${MPLAPACK_CUDA_ARCHITECTURES})
  list(REMOVE_AT _mplapack_mpc_cuda_extra_archs 0)
  foreach(_arch IN LISTS _mplapack_mpc_cuda_extra_archs)
    string(APPEND _mplapack_mpc_cuda_nvcc " -gencode arch=compute_${_arch},code=sm_${_arch}")
  endforeach()
  # mpc_cuda derives the CUDA library directory from NVCC, which now has flags
  set(_mplapack_mpc_cuda_libdir "${CUDAToolkit_LIBRARY_DIR}")
  # mpc_cuda builds from the GMP/MPFR/MPC sources of external/{gmp,mpfr,mpc}
  # and must not be built in parallel.
  find_program(MPLAPACK_PATCH_EXECUTABLE patch REQUIRED)
  set(_mplapack_mpc_cuda_patch "${_mplapack_external}/mpc_cuda/patches/patch-cu_freal-correctly-rounded-div-sqrt")
  ExternalProject_Add(mpc_cuda_bundled
    URL "${_mplapack_external}/mpc_cuda/download/mpc_cuda-0.0.1-abd9407.tar.xz"
    URL_HASH SHA256=bcac28655a62d381ffe38fc09704aca39f45610d1f697c351b4fc0acaa2278e7
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/external/mpc_cuda"
    PATCH_COMMAND
      ${CMAKE_COMMAND} -E tar xf "${_mplapack_external}/gmp/download/gmp-6.3.0.tar.xz"
      COMMAND ${CMAKE_COMMAND} -E tar xf "${_mplapack_external}/mpfr/download/mpfr-4.2.2.tar.bz2"
      COMMAND ${CMAKE_COMMAND} -E tar xf "${_mplapack_external}/mpc/download/mpc-1.4.1.tar.xz"
      COMMAND "${MPLAPACK_PATCH_EXECUTABLE}" -p1 --batch -i "${_mplapack_mpc_cuda_patch}"
    CONFIGURE_COMMAND <SOURCE_DIR>/configure "--prefix=${_mplapack_mpc_cuda_root}"
      "--with-cuda-arch=sm_${_mplapack_mpc_cuda_arch}" "NVCC=${_mplapack_mpc_cuda_nvcc}"
    BUILD_IN_SOURCE ON
    # C locale: tools/build_lib.sh prunes objects with sort | comm (see
    # external/mpc_cuda/Makefile.am)
    BUILD_COMMAND ${CMAKE_COMMAND} -E env MAKEFLAGS= LC_ALL=C "CUDA_LIBDIR=${_mplapack_mpc_cuda_libdir}" make -j1
    INSTALL_COMMAND ${CMAKE_COMMAND} -E env MAKEFLAGS= make -j1 install
    BUILD_BYPRODUCTS "${MPLAPACK_MPC_CUDA_LIBRARY}")
  # A changed patch invalidates download/extract, not just patch: never apply
  # an updated diff on top of an already patched source tree.
  ExternalProject_Add_Step(mpc_cuda_bundled patch_input
    INDEPENDENT TRUE
    DEPENDEES mkdir DEPENDERS download
    DEPENDS "${_mplapack_mpc_cuda_patch}")
endif()

if(MPLAPACK_ENABLE_MPFR_CUDA)
  # cuda/Rgemm.cpp shadows optimized/mpfr/Rgemm.cpp and the Level 1/2/3
  # cuda/R*.cpp the dispatchers of optimized/common/blas (whose CPU code they
  # call through common/blas/*_omp.cpp); the *_cuda files are helpers.
  set(_mplapack_mpfr_cuda_dir "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/mpfr/cuda")
  set(_mplapack_mpfr_cuda_l3 "")
  foreach(_r Rsyrk Rsyr2k Rsymm Rtrmm Rtrsm Rgemv Rsymv Rtrmv Rger Rsyr Rsyr2 Rscal Rswap Rrot)
    list(APPEND _mplapack_mpfr_cuda_l3 "${_mplapack_mpfr_cuda_dir}/${_r}.cpp")
  endforeach()
  # The complex BLAS (and Rgemmtr): the dispatchers of common/blas_cuda
  # (also those of dd/td/qd) on Cblas_mpfr_bridge_cuda.cpp.
  # The device code of the complex BLAS, split for parallel compilation;
  # MPLAPACK_MPFR_CUDA_RT_LARGE adds the 8192/16384-bit runtime-precision
  # types (their device code takes long to compile).
  set(_mplapack_mpfr_cuda_cblas_device "")
  foreach(_f Cblas_mpfr_device_cuda Cblas_mpfr_device_rt_cuda Cblas_mpfr_device_rt64_cuda
             Cblas_mpfr_device_rt128_cuda Cblas_mpfr_device_rt256_cuda)
    list(APPEND _mplapack_mpfr_cuda_cblas_device "${_mplapack_mpfr_cuda_dir}/${_f}.cu")
  endforeach()
  set(_mplapack_mpfr_cuda_cblas ${_mplapack_blas_cuda_sources})
  list(FILTER _mplapack_mpfr_cuda_cblas INCLUDE REGEX "/(C[A-Za-z0-9]+|Rgemmtr|Rsbmv|Rspmv|Rgbmv|Rtbmv|Rtpmv|Rspr|Rspr2)\\.cpp$")
  list(FILTER _mplapack_mpfr_cuda_cblas EXCLUDE REGEX "/Cblas_bridge_cuda\\.cpp$")
  mplapack_add_accelerator_backend(mpfr cuda ${_solver_mpfr_sources}
    "${_mplapack_mpfr_cuda_dir}/Rgetrf.cpp"
    "${_mplapack_mpfr_cuda_dir}/Cgetrf.cpp"
    "${_mplapack_mpfr_cuda_dir}/Rpotrf.cpp"
    "${_mplapack_mpfr_cuda_dir}/Cpotrf.cpp"
    "${_mplapack_mpfr_cuda_dir}/Cfactor_bridge_cuda.cpp"
    "${_mplapack_mpfr_cuda_dir}/Cfactor_device_cuda.cu"
    "${_mplapack_mpfr_cuda_dir}/Cfactor_rt_cuda.cu"
    "${_mplapack_mpfr_cuda_dir}/Rfactor_bridge_cuda.cpp"
    "${_mplapack_mpfr_cuda_dir}/Rfactor_device_cuda.cu"
    "${_mplapack_mpfr_cuda_dir}/Rfactor_rt_cuda.cu"
    "${_mplapack_mpfr_cuda_dir}/Rgemm.cpp"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_bridge_cuda.cpp"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_device_cuda.cu"
    "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu"
    ${_mplapack_mpfr_cuda_l3}
    "${_mplapack_mpfr_cuda_dir}/blas3_bridge_cuda.cpp"
    "${_mplapack_mpfr_cuda_dir}/blas3_device_cuda.cu"
    ${_mplapack_mpfr_cuda_cblas}
    "${_mplapack_mpfr_cuda_dir}/Cblas_mpfr_bridge_cuda.cpp"
    ${_mplapack_mpfr_cuda_cblas_device})
  if(MPLAPACK_MPFR_CUDA_RT_LARGE)
    target_compile_definitions(mplapack_mpfr_opt_cuda PRIVATE MPLAPACK_MPFR_CUDA_RT_LARGE)
  endif()
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

include(cmake/MplapackSolverTests.cmake)

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
  foreach(b IN LISTS MPLAPACK_BACKENDS)
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
  foreach(accelerator dd_opt_cuda qd_opt_cuda td_opt_cuda mpfr_opt_cuda binary128_opt_opencl)
    if(TARGET mplapack_${accelerator})
      add_test(NAME unique_symbols_${accelerator}
        COMMAND ${MPLAPACK_CHECK_UNIQUE_SYMBOLS_ENV}
                "${CMAKE_CURRENT_SOURCE_DIR}/misc/check_unique_symbols.sh"
                "$<TARGET_FILE:mplapack_${accelerator}>")
    endif()
  endforeach()
  # Shared by both test directories: executable-owned test overrides rely on the
  # linker keeping the first (overriding) definition of duplicated symbols.
  include(CheckLinkerFlag)
  check_linker_flag(CXX "-Wl,--allow-multiple-definition" MPLAPACK_LD_ALLOW_MULTIDEF)
  # The mpblas/mplapack correctness tests compare against a reference Fortran
  # BLAS/LAPACK and need its Fortran runtime.
  include(MplapackReferenceLapack)
  find_library(MPLAPACK_GFORTRAN_LIBRARY NAMES gfortran)
  add_subdirectory(mpblas/test)
  # Optimized CPU Rgemm of dd/td/qd/binary128: the blocked engine bit-identical
  # to the OpenMP loops; Winograd, Ozaki (and the branch-free td/qd kernels)
  # exact on integers and bounded on reals.
  set(_mplapack_opt_common "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/common")
  # Parallel Level 1/2/3 BLAS (real: Rblas_omp_test, complex: Cblas_omp_test)
  # of every optimized library (common/blas): bit-identical to the reference
  # routines.
  foreach(b IN LISTS MPLAPACK_BACKENDS)
    if(TARGET mplapack_${b}_opt)
      add_executable(${b}_opt_Rblas_omp_test "${_mplapack_opt_common}/test/Rblas_omp_test.cpp")
      target_link_libraries(${b}_opt_Rblas_omp_test PRIVATE mplapack_${b}_opt)
      add_test(NAME ${b}_opt_Rblas_omp COMMAND ${b}_opt_Rblas_omp_test)
      set_tests_properties(${b}_opt_Rblas_omp PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
      add_executable(${b}_opt_Cblas_omp_test "${_mplapack_opt_common}/test/Cblas_omp_test.cpp")
      target_link_libraries(${b}_opt_Cblas_omp_test PRIVATE mplapack_${b}_opt)
      add_test(NAME ${b}_opt_Cblas_omp COMMAND ${b}_opt_Cblas_omp_test)
      set_tests_properties(${b}_opt_Cblas_omp PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
      add_executable(${b}_opt_Rlapack_omp_test "${_mplapack_opt_common}/test/Rlapack_omp_test.cpp")
      target_link_libraries(${b}_opt_Rlapack_omp_test PRIVATE mplapack_${b}_opt)
      add_test(NAME ${b}_opt_Rlapack_omp COMMAND ${b}_opt_Rlapack_omp_test)
      set_tests_properties(${b}_opt_Rlapack_omp PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
      add_executable(${b}_opt_Rgemm_winograd_omp_test "${_mplapack_opt_common}/test/Rgemm_winograd_omp_test.cpp")
      target_link_libraries(${b}_opt_Rgemm_winograd_omp_test PRIVATE mplapack_${b}_opt)
      add_test(NAME ${b}_opt_Rgemm_winograd_omp COMMAND ${b}_opt_Rgemm_winograd_omp_test)
      set_tests_properties(${b}_opt_Rgemm_winograd_omp PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
    endif()
  endforeach()
  foreach(b dd td qd binary128)
    if(TARGET mplapack_${b}_opt)
      string(TOUPPER "${b}" _B)
      add_executable(${b}_opt_Rgemm_test "${_mplapack_opt_common}/test/Rgemm_opt_test.cpp")
      target_compile_definitions(${b}_opt_Rgemm_test PRIVATE RGEMM_TEST_${_B})
      target_link_libraries(${b}_opt_Rgemm_test PRIVATE mplapack_${b}_opt)
      add_test(NAME ${b}_opt_Rgemm COMMAND ${b}_opt_Rgemm_test)
    endif()
  endforeach()
  # Full factorization gates, CPU and GPU-BLAS-per-call variants.
  foreach(b IN LISTS MPLAPACK_BACKENDS)
    foreach(flavor opt opt_cuda)
      if(NOT TARGET mplapack_${b}_${flavor})
        continue()
      endif()
      set(_factor ${b}_${flavor}_Rfactor_test)
      add_executable(${_factor} mplapack/optimized/common/test/Rfactor_test.cpp)
      target_link_libraries(${_factor} PRIVATE mplapack_${b}_${flavor})
      mplapack_configure_dd_arithmetic(${_factor})
      if(flavor STREQUAL "opt_cuda")
        target_compile_definitions(${_factor} PRIVATE MPLAPACK_FACTOR_CUDA)
        target_link_libraries(${_factor} PRIVATE CUDA::cudart)
      endif()
      foreach(mode default fast)
        add_test(NAME ${b}_${flavor}_Rfactor_${mode} COMMAND ${_factor})
        if(mode STREQUAL "fast")
          set(_fast 1)
        else()
          set(_fast 0)
        endif()
        set_tests_properties(${b}_${flavor}_Rfactor_${mode} PROPERTIES
          ENVIRONMENT "OMP_NUM_THREADS=4;MPLAPACK_LAPACK_FAST=${_fast}"
          SKIP_RETURN_CODE 77 TIMEOUT 1800)
      endforeach()
    endforeach()
  endforeach()
  # The dd/td/qd arithmetic of the SIMD and GPU kernels against libQD3's
  # (branch-free) operators, bit for bit.
  foreach(b dd td qd)
    if(TARGET mplapack_${b}_opt)
      add_executable(qd_ops_test "${_mplapack_opt_common}/test/Rqd_ops_test.cpp")
      target_link_libraries(qd_ops_test PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(qd_ops_test)
      add_test(NAME qd_ops COMMAND qd_ops_test)
      break()
    endif()
  endforeach()
  # dd/td/qd CUDA Rgemm.  The host variant runs the GPU code on the host and
  # needs neither CUDA nor a GPU; the device variant is skipped (77) without
  # a device.
  foreach(b dd td qd)
    if(TARGET mplapack_${b}_opt)
      string(TOUPPER "${b}" _B)
      set(_exe ${b}_cuda_Rgetrf_host_test)
      add_executable(${_exe} mplapack/optimized/common/test/Rgetrf_gpu_test.cpp)
      target_compile_definitions(${_exe} PRIVATE MPLAPACK_GEMM_CUDA_HOST_EMULATION)
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(${_exe})
      add_test(NAME ${b}_cuda_Rgetrf_host COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rgetrf_host PROPERTIES TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=1")
      set(_exe ${b}_cuda_Rpotrf_host_test)
      add_executable(${_exe} mplapack/optimized/common/test/Rpotrf_gpu_test.cpp)
      target_compile_definitions(${_exe} PRIVATE MPLAPACK_GEMM_CUDA_HOST_EMULATION)
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(${_exe})
      add_test(NAME ${b}_cuda_Rpotrf_host COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rpotrf_host PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
      set(_exe ${b}_cuda_Rlapack_host_test)
      add_executable(${_exe} "${_mplapack_opt_common}/test/Rlapack_gpu_test.cpp")
      target_compile_definitions(${_exe} PRIVATE MPLAPACK_GEMM_CUDA_HOST_EMULATION)
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(${_exe})
      add_test(NAME ${b}_cuda_Rlapack_host COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rlapack_host PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
      set(_exe ${b}_cuda_Rgemm_host_test)
      add_executable(${_exe}
        "${_mplapack_opt_common}/test/Rgemm_gpu_test.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${b}/cuda/Rgemm_bridge_cuda.cpp")
      target_compile_definitions(${_exe} PRIVATE RGEMM_TEST_${_B} MPLAPACK_GEMM_CUDA_HOST_EMULATION)
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(${_exe})
      add_test(NAME ${b}_cuda_Rgemm_host COMMAND ${_exe})
      # Level 1/2/3 BLAS of common/blas_cuda, every call on the (emulated) GPU
      set(_exe ${b}_cuda_Rblas_host_test)
      add_executable(${_exe} "${_mplapack_opt_common}/test/Rblas_gpu_test.cpp" ${_mplapack_blas_cuda_sources})
      target_compile_definitions(${_exe} PRIVATE MPLAPACK_GEMM_CUDA_HOST_EMULATION)
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(${_exe})
      add_test(NAME ${b}_cuda_Rblas_host COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rblas_host PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
      set(_exe ${b}_cuda_Cblas_host_test)
      add_executable(${_exe} "${_mplapack_opt_common}/test/Cblas_gpu_test.cpp" ${_mplapack_blas_cuda_sources})
      target_compile_definitions(${_exe} PRIVATE MPLAPACK_GEMM_CUDA_HOST_EMULATION)
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt)
      mplapack_configure_dd_arithmetic(${_exe})
      add_test(NAME ${b}_cuda_Cblas_host COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Cblas_host PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
    endif()
    if(TARGET mplapack_${b}_opt_cuda)
      set(_exe ${b}_cuda_Rgetrf_device_test)
      add_executable(${_exe} mplapack/optimized/common/test/Rgetrf_gpu_device_test.cu)
      target_compile_options(${_exe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>" "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-ffp-contract=off>")
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt_cuda CUDA::cudart)
      mplapack_configure_dd_arithmetic(${_exe})
      set_target_properties(${_exe} PROPERTIES CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}" CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
      add_test(NAME ${b}_cuda_Rgetrf_device COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rgetrf_device PROPERTIES TIMEOUT 1800 SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1;MPLAPACK_${_B}_CUDA=0")
      set(_exe ${b}_cuda_Rpotrf_device_test)
      add_executable(${_exe} mplapack/optimized/common/test/Rpotrf_gpu_device_test.cu)
      target_compile_options(${_exe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>" "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-ffp-contract=off>")
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt_cuda CUDA::cudart)
      mplapack_configure_dd_arithmetic(${_exe})
      set_target_properties(${_exe} PROPERTIES CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}" CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
      add_test(NAME ${b}_cuda_Rpotrf_device COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rpotrf_device PROPERTIES SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1;MPLAPACK_${_B}_CUDA=0")
      set(_exe ${b}_cuda_Rlapack_device_test)
      add_executable(${_exe} "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${b}/cuda/Rlapack_device_test.cu")
      # This CUDA TU also contains the host libQD3 oracle.
      target_compile_options(${_exe} PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>"
        "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-ffp-contract=off>")
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt CUDA::cudart CUDA::cublas)
      mplapack_configure_dd_arithmetic(${_exe})
      set_target_properties(${_exe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
      add_test(NAME ${b}_cuda_Rlapack_device COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rlapack_device PROPERTIES SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1")
      set(_exe ${b}_cuda_Rgemm_device_test)
      add_executable(${_exe}
        "${_mplapack_opt_common}/test/Rgemm_gpu_test.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${b}/cuda/Rgemm_bridge_cuda.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/mpblas/optimized/${b}/cuda/Rgemm_gpu_cuda.cu")
      target_compile_definitions(${_exe} PRIVATE RGEMM_TEST_${_B})
      target_compile_options(${_exe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt CUDA::cudart CUDA::cublas)
      mplapack_configure_dd_arithmetic(${_exe})
      set_target_properties(${_exe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON)
      add_test(NAME ${b}_cuda_Rgemm_device COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rgemm_device PROPERTIES SKIP_RETURN_CODE 77)
      set(_exe ${b}_cuda_Rblas_device_test)
      add_executable(${_exe} "${_mplapack_opt_common}/test/Rblas_gpu_test.cpp")
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt_cuda CUDA::cudart)
      add_test(NAME ${b}_cuda_Rblas_device COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Rblas_device PROPERTIES SKIP_RETURN_CODE 77)
      set(_exe ${b}_cuda_Cblas_device_test)
      add_executable(${_exe} "${_mplapack_opt_common}/test/Cblas_gpu_test.cpp")
      target_link_libraries(${_exe} PRIVATE mplapack_${b}_opt_cuda CUDA::cudart)
      add_test(NAME ${b}_cuda_Cblas_device COMMAND ${_exe})
      set_tests_properties(${b}_cuda_Cblas_device PROPERTIES SKIP_RETURN_CODE 77)
    endif()
  endforeach()
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
    add_executable(mpfr_cpu_Rfactor_bench "${_mplapack_mpfr_cuda_dir}/test/Rfactor_bench.cpp")
    target_compile_definitions(mpfr_cpu_Rfactor_bench PRIVATE MPLAPACK_FACTOR_BENCH_CPU_ONLY)
    target_link_libraries(mpfr_cpu_Rfactor_bench PRIVATE mplapack_mpfr_opt)
    add_executable(mpfr_cuda_Rfactor_bench "${_mplapack_mpfr_cuda_dir}/test/Rfactor_bench.cpp")
    target_link_libraries(mpfr_cuda_Rfactor_bench PRIVATE mplapack_mpfr_opt_cuda CUDA::cudart)
    add_executable(mpfr_cpu_Cfactor_bench "${_mplapack_mpfr_cuda_dir}/test/Cfactor_bench.cpp")
    target_compile_definitions(mpfr_cpu_Cfactor_bench PRIVATE MPLAPACK_FACTOR_BENCH_CPU_ONLY)
    target_link_libraries(mpfr_cpu_Cfactor_bench PRIVATE mplapack_mpfr_opt)
    add_executable(mpfr_cuda_Cfactor_bench "${_mplapack_mpfr_cuda_dir}/test/Cfactor_bench.cpp")
    target_link_libraries(mpfr_cuda_Cfactor_bench PRIVATE mplapack_mpfr_opt_cuda CUDA::cudart)
    add_executable(mpfr_cuda_Rfactor_panel_bench
      "${_mplapack_mpfr_cuda_dir}/test/Rfactor_panel_bench.cpp"
      "${_mplapack_mpfr_cuda_dir}/Rfactor_device_cuda.cu"
      "${_mplapack_mpfr_cuda_dir}/Rfactor_rt_cuda.cu")
    target_include_directories(mpfr_cuda_Rfactor_panel_bench PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
    target_link_libraries(mpfr_cuda_Rfactor_panel_bench PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
    target_compile_options(mpfr_cuda_Rfactor_panel_bench PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
    set_target_properties(mpfr_cuda_Rfactor_panel_bench PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
      CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON
      CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
    if(_mplapack_mpc_cuda_bundled)
      add_dependencies(mpfr_cuda_Rfactor_panel_bench mpc_cuda_bundled)
    endif()
    add_executable(mpfr_cuda_Cfactor_panel_bench
      "${_mplapack_mpfr_cuda_dir}/test/Cfactor_panel_bench.cpp"
      "${_mplapack_mpfr_cuda_dir}/Cfactor_device_cuda.cu"
      "${_mplapack_mpfr_cuda_dir}/Cfactor_rt_cuda.cu"
      "${_mplapack_mpfr_cuda_dir}/test/Cfactor_device_globals.cu")
    target_include_directories(mpfr_cuda_Cfactor_panel_bench PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
    target_link_libraries(mpfr_cuda_Cfactor_panel_bench PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
    target_compile_options(mpfr_cuda_Cfactor_panel_bench PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
    set_target_properties(mpfr_cuda_Cfactor_panel_bench PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
      CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON
      CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
    if(_mplapack_mpc_cuda_bundled)
      add_dependencies(mpfr_cuda_Cfactor_panel_bench mpc_cuda_bundled)
    endif()
    if(UNIX AND NOT APPLE)
      foreach(_factor_bench mpfr_cpu_Rfactor_bench mpfr_cuda_Rfactor_bench mpfr_cuda_Rfactor_panel_bench mpfr_cpu_Cfactor_bench mpfr_cuda_Cfactor_bench mpfr_cuda_Cfactor_panel_bench)
        target_link_options(${_factor_bench} PRIVATE "-Wl,--disable-new-dtags")
      endforeach()
    endif()
    add_executable(mpfr_cuda_Rtrsm_round_bench "${_mplapack_mpfr_cuda_dir}/test/Rtrsm_round_bench.cpp")
    target_link_libraries(mpfr_cuda_Rtrsm_round_bench PRIVATE mplapack_mpfr_opt_cuda CUDA::cudart)
    add_executable(mpfr_cuda_cu_fround_bench "${_mplapack_mpfr_cuda_dir}/test/cu_fround_bench.cu")
    target_include_directories(mpfr_cuda_cu_fround_bench PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
    target_link_libraries(mpfr_cuda_cu_fround_bench PRIVATE "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
    target_compile_options(mpfr_cuda_cu_fround_bench PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
    set_target_properties(mpfr_cuda_cu_fround_bench PROPERTIES
      CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON CUDA_SEPARABLE_COMPILATION ON
      CUDA_RESOLVE_DEVICE_SYMBOLS ON CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
    if(_mplapack_mpc_cuda_bundled)
      add_dependencies(mpfr_cuda_cu_fround_bench mpc_cuda_bundled)
    endif()
    foreach(_variant host device)
      set(_fexe mpfr_cuda_Rfactor_${_variant}_test)
      if(_variant STREQUAL "host")
        set(_factor_engine_sources
          "${_mplapack_mpfr_cuda_dir}/test/Rfactor_fixed_host.cpp"
          "${_mplapack_mpfr_cuda_dir}/test/Rfactor_runtime_host.cpp")
      else()
        set(_factor_engine_sources
          "${_mplapack_mpfr_cuda_dir}/Rfactor_device_cuda.cu"
          "${_mplapack_mpfr_cuda_dir}/Rfactor_rt_cuda.cu")
      endif()
      add_executable(${_fexe}
        "${_mplapack_mpfr_cuda_dir}/test/Rfactor_cuda_test.cpp"
        "${_mplapack_mpfr_cuda_dir}/Rfactor_bridge_cuda.cpp"
        ${_factor_engine_sources})
      target_include_directories(${_fexe} PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
      target_link_libraries(${_fexe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      if(UNIX AND NOT APPLE)
        target_link_options(${_fexe} PRIVATE "-Wl,--disable-new-dtags")
      endif()
      if(_variant STREQUAL "host")
        target_compile_definitions(${_fexe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
      else()
        target_link_libraries(${_fexe} PRIVATE "${MPLAPACK_MPC_CUDA_LIBRARY}")
        target_compile_options(${_fexe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
        set_target_properties(${_fexe} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
          CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON
          CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
      endif()
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_fexe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_Rfactor_${_variant} COMMAND ${_fexe})
      add_test(NAME mpfr_cuda_Rfactor_${_variant}_runtime COMMAND ${_fexe})
      set_tests_properties(mpfr_cuda_Rfactor_${_variant} PROPERTIES SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1")
      set_tests_properties(mpfr_cuda_Rfactor_${_variant}_runtime PROPERTIES
        SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1;MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1")
      set(_cfexe mpfr_cuda_Cfactor_${_variant}_test)
      if(_variant STREQUAL "host")
        set(_cfactor_sources
          "${_mplapack_mpfr_cuda_dir}/test/Cfactor_fixed_host.cpp"
          "${_mplapack_mpfr_cuda_dir}/test/Cfactor_runtime_host.cpp")
      else()
        set(_cfactor_sources
          "${_mplapack_mpfr_cuda_dir}/Cfactor_device_cuda.cu"
          "${_mplapack_mpfr_cuda_dir}/Cfactor_rt_cuda.cu"
          "${_mplapack_mpfr_cuda_dir}/test/Cfactor_device_globals.cu")
      endif()
      add_executable(${_cfexe}
        "${_mplapack_mpfr_cuda_dir}/test/Cfactor_cuda_test.cpp"
        "${_mplapack_mpfr_cuda_dir}/Cfactor_bridge_cuda.cpp" ${_cfactor_sources})
      target_include_directories(${_cfexe} PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
      target_link_libraries(${_cfexe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      if(UNIX AND NOT APPLE)
        target_link_options(${_cfexe} PRIVATE "-Wl,--disable-new-dtags")
      endif()
      if(_variant STREQUAL "host")
        target_compile_definitions(${_cfexe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
      else()
        target_link_libraries(${_cfexe} PRIVATE "${MPLAPACK_MPC_CUDA_LIBRARY}")
        target_compile_options(${_cfexe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
        set_target_properties(${_cfexe} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
          CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON
          CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
      endif()
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_cfexe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_Cfactor_${_variant} COMMAND ${_cfexe})
      add_test(NAME mpfr_cuda_Cfactor_${_variant}_runtime COMMAND ${_cfexe})
      set_tests_properties(mpfr_cuda_Cfactor_${_variant} PROPERTIES SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1")
      set_tests_properties(mpfr_cuda_Cfactor_${_variant}_runtime PROPERTIES SKIP_RETURN_CODE 77
        ENVIRONMENT "OMP_NUM_THREADS=1;MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1")
      set(_roundexe mpfr_cuda_cu_fround_${_variant}_test)
      if(_variant STREQUAL "host")
        add_executable(${_roundexe} "${_mplapack_mpfr_cuda_dir}/test/cu_fround_test.cpp")
      else()
        add_executable(${_roundexe} "${_mplapack_mpfr_cuda_dir}/test/cu_fround_device_test.cu")
        target_link_libraries(${_roundexe} PRIVATE CUDA::cudart)
        target_compile_options(${_roundexe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
        set_target_properties(${_roundexe} PROPERTIES CUDA_STANDARD 17
          CUDA_STANDARD_REQUIRED ON CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}")
      endif()
      target_include_directories(${_roundexe} PRIVATE "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}")
      target_link_libraries(${_roundexe} PRIVATE mplapack_mpfr_opt)
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_roundexe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_cu_fround_${_variant} COMMAND ${_roundexe})
      set_tests_properties(mpfr_cuda_cu_fround_${_variant} PROPERTIES SKIP_RETURN_CODE 77)
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

      # Level 3 (Rsyrk, Rsyr2k, Rsymm, Rtrmm, Rtrsm): bit for bit against the
      # reference implementation, on both engines.
      set(_l3exe mpfr_cuda_blas3_${_variant}_test)
      add_executable(${_l3exe}
        "${_mplapack_mpfr_cuda_dir}/test/blas3_cuda_test.cpp"
        "${_mplapack_mpfr_cuda_dir}/blas3_bridge_cuda.cpp"
        "${_mplapack_mpfr_cuda_dir}/blas3_device_cuda.cu"
        "${_fixed}"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu")
      if(_variant STREQUAL "host")
        target_compile_definitions(${_l3exe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
        target_link_libraries(${_l3exe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      else()
        target_link_libraries(${_l3exe} PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
      endif()
      target_include_directories(${_l3exe} PRIVATE
        "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}" "${_mplapack_mpfr_cuda_dir}")
      target_compile_options(${_l3exe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
      set_target_properties(${_l3exe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON)
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_l3exe} mpc_cuda_bundled)
      endif()
      # Level 1/2/3 through the dispatchers, every call on the GPU
      # (common/test/Rblas_gpu_test.cpp, cu_mpfr engine).
      set(_bexe mpfr_cuda_Rblas_${_variant}_test)
      add_executable(${_bexe}
        "${_mplapack_opt_common}/test/Rblas_gpu_test.cpp"
        ${_mplapack_mpfr_cuda_l3}
        "${_mplapack_mpfr_cuda_dir}/blas3_bridge_cuda.cpp"
        "${_mplapack_mpfr_cuda_dir}/blas3_device_cuda.cu"
        "${_fixed}"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu")
      if(_variant STREQUAL "host")
        target_compile_definitions(${_bexe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
        target_link_libraries(${_bexe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      else()
        target_link_libraries(${_bexe} PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
      endif()
      target_include_directories(${_bexe} PRIVATE
        "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}" "${_mplapack_mpfr_cuda_dir}")
      target_compile_options(${_bexe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
      set_target_properties(${_bexe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON)
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_bexe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_Rblas_${_variant} COMMAND ${_bexe})
      set_tests_properties(mpfr_cuda_Rblas_${_variant} PROPERTIES SKIP_RETURN_CODE 77)

      # Complex Level 1/2/3 (and Rgemmtr) through the dispatchers, every
      # call on the GPU (common/test/Cblas_gpu_test.cpp, cu_freal engine).
      set(_cexe mpfr_cuda_Cblas_${_variant}_test)
      add_executable(${_cexe}
        "${_mplapack_opt_common}/test/Cblas_gpu_test.cpp"
        ${_mplapack_mpfr_cuda_cblas}
        "${_mplapack_mpfr_cuda_dir}/Cblas_mpfr_bridge_cuda.cpp"
        ${_mplapack_mpfr_cuda_cblas_device}
        "${_fixed}"
        "${_mplapack_mpfr_cuda_dir}/Rgemm_rt_cuda.cu")
      if(_variant STREQUAL "host")
        target_compile_definitions(${_cexe} PRIVATE MPLAPACK_MPFR_CUDA_HOST_EMULATION)
        target_link_libraries(${_cexe} PRIVATE mplapack_mpfr_opt CUDA::cudart)
      else()
        target_link_libraries(${_cexe} PRIVATE mplapack_mpfr_opt "${MPLAPACK_MPC_CUDA_LIBRARY}" CUDA::cudart)
      endif()
      if(MPLAPACK_MPFR_CUDA_RT_LARGE)
        target_compile_definitions(${_cexe} PRIVATE MPLAPACK_MPFR_CUDA_RT_LARGE)
      endif()
      target_include_directories(${_cexe} PRIVATE
        "${MPLAPACK_MPC_CUDA_INCLUDE_DIR}" "${_mplapack_mpfr_cuda_dir}")
      target_compile_options(${_cexe} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
      set_target_properties(${_cexe} PROPERTIES
        CUDA_ARCHITECTURES "${MPLAPACK_CUDA_ARCHITECTURES}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON)
      if(_mplapack_mpc_cuda_bundled)
        add_dependencies(${_cexe} mpc_cuda_bundled)
      endif()
      add_test(NAME mpfr_cuda_Cblas_${_variant} COMMAND ${_cexe})
      set_tests_properties(mpfr_cuda_Cblas_${_variant} PROPERTIES SKIP_RETURN_CODE 77)
      # the cu_freal engine (the test forces the runtime one otherwise)
      add_test(NAME mpfr_cuda_Cblas_${_variant}_fixed COMMAND ${_cexe})
      set_tests_properties(mpfr_cuda_Cblas_${_variant}_fixed PROPERTIES
        SKIP_RETURN_CODE 77 ENVIRONMENT "MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0")

      add_test(NAME mpfr_cuda_blas3_${_variant} COMMAND ${_l3exe})
      add_test(NAME mpfr_cuda_blas3_${_variant}_runtime COMMAND ${_l3exe})
      set_tests_properties(mpfr_cuda_blas3_${_variant} PROPERTIES SKIP_RETURN_CODE 77)
      set_tests_properties(mpfr_cuda_blas3_${_variant}_runtime PROPERTIES
        SKIP_RETURN_CODE 77 ENVIRONMENT "MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1")
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
