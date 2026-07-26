# FindGMP.cmake — locate the GNU MP C interface.
#
# Provides GMP::GMP and the classic variables GMP_FOUND,
# GMP_INCLUDE_DIRS, GMP_LIBRARIES, and GMP_VERSION.

find_package(PkgConfig QUIET)
pkg_check_modules(PC_GMP QUIET gmp)

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  HINTS ${PC_GMP_INCLUDEDIR} ${PC_GMP_INCLUDE_DIRS})

find_library(GMP_LIBRARY
  NAMES gmp
  HINTS ${PC_GMP_LIBDIR} ${PC_GMP_LIBRARY_DIRS})

set(GMP_VERSION ${PC_GMP_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR
  VERSION_VAR GMP_VERSION)

if(GMP_FOUND)
  set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
  set(GMP_LIBRARIES ${GMP_LIBRARY})
  if(NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(GMP::GMP PROPERTIES
      IMPORTED_LOCATION "${GMP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
