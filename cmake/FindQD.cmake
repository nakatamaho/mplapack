# FindQD.cmake — locate the QD library (double-double / quad-double).
#
# QD ships a qd_real.h under a qd/ subdirectory and a libqd. Some installs
# also provide a qd.pc; we use it as a hint when present.
#
# Provides imported target QD::QD and variables
# QD_FOUND, QD_INCLUDE_DIRS, QD_LIBRARIES, QD_VERSION.

find_package(PkgConfig QUIET)
pkg_check_modules(PC_QD QUIET qd)

# qd_real.h is included as <qd/qd_real.h>, so the include root is the directory
# that *contains* the qd/ subdirectory.
find_path(QD_INCLUDE_DIR
  NAMES qd/qd_real.h
  HINTS ${PC_QD_INCLUDEDIR} ${PC_QD_INCLUDE_DIRS})

find_library(QD_LIBRARY
  NAMES qd
  HINTS ${PC_QD_LIBDIR} ${PC_QD_LIBRARY_DIRS})

set(QD_VERSION ${PC_QD_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QD
  REQUIRED_VARS QD_LIBRARY QD_INCLUDE_DIR
  VERSION_VAR QD_VERSION)

if(QD_FOUND)
  set(QD_INCLUDE_DIRS ${QD_INCLUDE_DIR})
  set(QD_LIBRARIES ${QD_LIBRARY})
  if(NOT TARGET QD::QD)
    add_library(QD::QD UNKNOWN IMPORTED)
    set_target_properties(QD::QD PROPERTIES
      IMPORTED_LOCATION "${QD_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${QD_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(QD_INCLUDE_DIR QD_LIBRARY)
