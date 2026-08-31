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

# Do not use a host qd.pc for a different custom QD library selected by the
# caller.  If the selected library has no matching pkg-config metadata, the
# fallback below restores the dependency information that can be inferred
# from the public QD archive.
set(QD_PKGCONFIG_FOUND FALSE)
if(PC_QD_FOUND)
  get_filename_component(_qd_selected_libdir "${QD_LIBRARY}" DIRECTORY)
  set(_qd_pc_libdirs ${PC_QD_LIBDIR} ${PC_QD_LIBRARY_DIRS}
      ${PC_QD_STATIC_LIBRARY_DIRS})
  list(REMOVE_DUPLICATES _qd_pc_libdirs)
  if(NOT _qd_pc_libdirs)
    set(QD_PKGCONFIG_FOUND TRUE)
  else()
    foreach(_qd_pc_libdir IN LISTS _qd_pc_libdirs)
      get_filename_component(_qd_pc_libdir_abs "${_qd_pc_libdir}" REALPATH)
      get_filename_component(_qd_selected_libdir_abs "${_qd_selected_libdir}" REALPATH)
      if(_qd_pc_libdir_abs STREQUAL _qd_selected_libdir_abs)
        set(QD_PKGCONFIG_FOUND TRUE)
      endif()
    endforeach()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QD
  REQUIRED_VARS QD_LIBRARY QD_INCLUDE_DIR
  VERSION_VAR QD_VERSION)

if(QD_FOUND)
  set(QD_INCLUDE_DIRS ${QD_INCLUDE_DIR})
  if(QD_PKGCONFIG_FOUND)
    list(APPEND QD_INCLUDE_DIRS ${PC_QD_INCLUDE_DIRS})
  endif()
  set(QD_LIBRARIES ${QD_LIBRARY})
  if(NOT TARGET QD::QD)
    add_library(QD::QD UNKNOWN IMPORTED)
    if(QD_PKGCONFIG_FOUND)
      set(_qd_interface_includes ${QD_INCLUDE_DIR} ${PC_QD_INCLUDE_DIRS}
          ${PC_QD_STATIC_INCLUDE_DIRS})
      set(_qd_interface_libs ${PC_QD_LIBRARIES}
          ${PC_QD_STATIC_LIBRARIES})
      set(_qd_interface_link_dirs ${PC_QD_LIBRARY_DIRS}
          ${PC_QD_STATIC_LIBRARY_DIRS})
      set(_qd_interface_compile_options ${PC_QD_CFLAGS_OTHER}
          ${PC_QD_STATIC_CFLAGS_OTHER})
      set(_qd_interface_link_options ${PC_QD_LDFLAGS_OTHER}
          ${PC_QD_STATIC_LDFLAGS_OTHER})
    else()
      set(_qd_interface_includes ${QD_INCLUDE_DIR})
      set(_qd_interface_libs)
      set(_qd_interface_link_dirs)
      set(_qd_interface_compile_options)
      set(_qd_interface_link_options)
    endif()
    list(REMOVE_ITEM _qd_interface_libs qd)
    list(REMOVE_DUPLICATES _qd_interface_libs)
    if(NOT QD_PKGCONFIG_FOUND AND UNIX)
      # QD's static archive uses libm, but a system without qd.pc gives us no
      # pkg-config closure from which to recover that dependency.
      list(APPEND _qd_interface_libs m)
    endif()
    set_target_properties(QD::QD PROPERTIES
      IMPORTED_LOCATION "${QD_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${_qd_interface_includes}"
      INTERFACE_LINK_LIBRARIES "${_qd_interface_libs}"
      INTERFACE_LINK_DIRECTORIES "${_qd_interface_link_dirs}"
      INTERFACE_COMPILE_OPTIONS "${_qd_interface_compile_options}"
      INTERFACE_LINK_OPTIONS "${_qd_interface_link_options}")
    list(APPEND QD_INCLUDE_DIRS ${_qd_interface_includes})
    list(APPEND QD_LIBRARIES ${_qd_interface_libs})
  endif()
endif()

mark_as_advanced(QD_INCLUDE_DIR QD_LIBRARY)
