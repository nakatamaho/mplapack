# FindMPFR.cmake — locate the MPFR library.
#
# Provides imported target MPFR::MPFR (links GMP::GMP) and variables
# MPFR_FOUND, MPFR_INCLUDE_DIRS, MPFR_LIBRARIES, MPFR_VERSION.

find_package(PkgConfig QUIET)
pkg_check_modules(PC_MPFR QUIET mpfr)

find_path(MPFR_INCLUDE_DIR
  NAMES mpfr.h
  HINTS ${PC_MPFR_INCLUDEDIR} ${PC_MPFR_INCLUDE_DIRS})

find_library(MPFR_LIBRARY
  NAMES mpfr
  HINTS ${PC_MPFR_LIBDIR} ${PC_MPFR_LIBRARY_DIRS})

set(MPFR_VERSION ${PC_MPFR_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
  REQUIRED_VARS MPFR_LIBRARY MPFR_INCLUDE_DIR
  VERSION_VAR MPFR_VERSION)

if(MPFR_FOUND)
  set(MPFR_INCLUDE_DIRS ${MPFR_INCLUDE_DIR})
  set(MPFR_LIBRARIES ${MPFR_LIBRARY})
  find_package(GMP QUIET)
  if(NOT TARGET MPFR::MPFR)
    add_library(MPFR::MPFR UNKNOWN IMPORTED)
    set_target_properties(MPFR::MPFR PROPERTIES
      IMPORTED_LOCATION "${MPFR_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}")
    if(TARGET GMP::GMP)
      set_target_properties(MPFR::MPFR PROPERTIES
        INTERFACE_LINK_LIBRARIES GMP::GMP)
    endif()
  endif()
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)
