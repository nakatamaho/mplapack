# Reference BLAS/LAPACK for the correctness tests.
#
# The mpblas and compare tests check MPLAPACK against the reference Fortran
# BLAS/LAPACK it was translated from.  An optimized BLAS is not a valid
# oracle: OpenBLAS, for example, returns 0 from dnrm2 for incx = 0 where
# reference BLAS 3.12 returns sqrt(n)*|x(1)|.  As external/lapack does for
# Autotools, build the bundled reference LAPACK tarball with the same
# patches, and link the tests against it.
#
# Sets BLAS_FOUND, LAPACK_FOUND, BLAS_LIBRARIES and LAPACK_LIBRARIES.  With
# MPLAPACK_TEST_BUNDLED_LAPACK=OFF, or without a Fortran compiler or the
# patch program, falls back to find_package(BLAS/LAPACK).

option(MPLAPACK_TEST_BUNDLED_LAPACK
  "Build the bundled reference LAPACK as the test oracle" ON)

set(_mplapack_ref_lapack_ok FALSE)
if(MPLAPACK_TEST_BUNDLED_LAPACK)
  include(CheckLanguage)
  check_language(Fortran)
  find_program(MPLAPACK_PATCH_EXECUTABLE patch)
  if(NOT CMAKE_Fortran_COMPILER)
    message(STATUS "Reference LAPACK: no Fortran compiler, using find_package(BLAS/LAPACK)")
  elseif(NOT MPLAPACK_PATCH_EXECUTABLE)
    message(STATUS "Reference LAPACK: no patch program, using find_package(BLAS/LAPACK)")
  else()
    set(_mplapack_ref_lapack_ok TRUE)
  endif()
endif()

if(NOT _mplapack_ref_lapack_ok)
  find_package(BLAS)
  find_package(LAPACK)
  return()
endif()

# Version, tarball checksum and patch list come from the Autotools files, so
# both build systems test against the same reference.
set(_ext "${PROJECT_SOURCE_DIR}/external")
file(READ "${_ext}/lapack/Makefile.am" _am)
if(NOT _am MATCHES "LAPACK_VERSION[ \t]*:=[ \t]*\\$\\(or \\$\\(LAPACK_VERSION\\),([0-9.]+)\\)")
  message(FATAL_ERROR "Cannot find the LAPACK version in external/lapack/Makefile.am")
endif()
set(_ver "${CMAKE_MATCH_1}")
string(REPLACE "." "_" _key "${_ver}")
if(NOT _am MATCHES "PATCHES_${_key} =[ \t]*\\\\\n(([ \t]+[^\n]*\n)+)")
  message(FATAL_ERROR "Cannot find PATCHES_${_key} in external/lapack/Makefile.am")
endif()
string(REGEX REPLACE "[ \t\\\\\n]+" ";" _patches "${CMAKE_MATCH_1}")
list(REMOVE_ITEM _patches "")
set(_patch_dir "${_ext}/lapack/patches/${_ver}")
# Same extra patches, in the same order, as the Autotools rule.
set(_extra "")
if(EXISTS "${_patch_dir}/patch-print")
  list(APPEND _extra patch-print)
endif()
if(MINGW AND EXISTS "${_patch_dir}/patch-mingw")
  list(APPEND _extra patch-mingw)
endif()
list(PREPEND _patches ${_extra})
list(JOIN _patches "|" _patches_arg)

set(_tarball "lapack/download/lapack-${_ver}.tar.gz")
file(STRINGS "${_ext}/distfiles.sha256" _sha_line REGEX " ${_tarball}$")
if(NOT _sha_line MATCHES "^([0-9a-f]+) ")
  message(FATAL_ERROR "No SHA256 for ${_tarball} in external/distfiles.sha256")
endif()
set(_sha "${CMAKE_MATCH_1}")

include(ExternalProject)
set(_prefix "${CMAKE_BINARY_DIR}/reference-lapack")
set(_blas "${_prefix}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(_lapack "${_prefix}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}")
ExternalProject_Add(mplapack_reference_lapack
  URL "${_ext}/${_tarball}"
  URL_HASH "SHA256=${_sha}"
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  PATCH_COMMAND "${CMAKE_COMMAND}"
    "-DSOURCE_DIR=<SOURCE_DIR>" "-DPATCH_DIR=${_patch_dir}"
    "-DPATCHES=${_patches_arg}" "-DPATCH=${MPLAPACK_PATCH_EXECUTABLE}"
    -P "${PROJECT_SOURCE_DIR}/cmake/ApplyLapackPatches.cmake"
  CMAKE_ARGS
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
    -DCMAKE_BUILD_TYPE=Release
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DBUILD_TESTING=OFF
    -DCBLAS=OFF
    -DLAPACKE=OFF
    "-DCMAKE_INSTALL_PREFIX=${_prefix}"
    -DCMAKE_INSTALL_LIBDIR=lib
  BUILD_BYPRODUCTS "${_blas}" "${_lapack}")

add_library(mplapack_reference_blas STATIC IMPORTED GLOBAL)
set_target_properties(mplapack_reference_blas PROPERTIES IMPORTED_LOCATION "${_blas}")
add_dependencies(mplapack_reference_blas mplapack_reference_lapack)
add_library(mplapack_reference_lapack_lib STATIC IMPORTED GLOBAL)
set_target_properties(mplapack_reference_lapack_lib PROPERTIES IMPORTED_LOCATION "${_lapack}")
add_dependencies(mplapack_reference_lapack_lib mplapack_reference_lapack)

set(BLAS_FOUND TRUE)
set(LAPACK_FOUND TRUE)
set(BLAS_LIBRARIES mplapack_reference_blas)
set(LAPACK_LIBRARIES mplapack_reference_lapack_lib)
message(STATUS "Reference LAPACK: bundled ${_ver} with ${_patches_arg}")
