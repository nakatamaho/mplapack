# FindMPC.cmake — locate the MPC library (complex arithmetic atop MPFR).
#
# Provides imported target MPC::MPC (links MPFR::MPFR) and variables
# MPC_FOUND, MPC_INCLUDE_DIRS, MPC_LIBRARIES, MPC_VERSION.

find_package(PkgConfig QUIET)
pkg_check_modules(PC_MPC QUIET mpc)

find_path(MPC_INCLUDE_DIR
  NAMES mpc.h
  HINTS ${PC_MPC_INCLUDEDIR} ${PC_MPC_INCLUDE_DIRS})

find_library(MPC_LIBRARY
  NAMES mpc
  HINTS ${PC_MPC_LIBDIR} ${PC_MPC_LIBRARY_DIRS})

set(MPC_VERSION ${PC_MPC_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPC
  REQUIRED_VARS MPC_LIBRARY MPC_INCLUDE_DIR
  VERSION_VAR MPC_VERSION)

if(MPC_FOUND)
  set(MPC_INCLUDE_DIRS ${MPC_INCLUDE_DIR})
  set(MPC_LIBRARIES ${MPC_LIBRARY})
  find_package(MPFR QUIET)
  if(NOT TARGET MPC::MPC)
    add_library(MPC::MPC UNKNOWN IMPORTED)
    set_target_properties(MPC::MPC PROPERTIES
      IMPORTED_LOCATION "${MPC_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MPC_INCLUDE_DIR}")
    if(TARGET MPFR::MPFR)
      set_target_properties(MPC::MPC PROPERTIES
        INTERFACE_LINK_LIBRARIES MPFR::MPFR)
    endif()
  endif()
endif()

mark_as_advanced(MPC_INCLUDE_DIR MPC_LIBRARY)
