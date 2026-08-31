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

# A globally discoverable gmp.pc is usable only when it describes the GMP
# library selected below.  This matters when a caller supplies a custom GMP
# prefix that has no pkg-config file but the host has another gmp.pc.
set(GMP_PKGCONFIG_FOUND FALSE)
if(PC_GMP_FOUND)
  get_filename_component(_gmp_selected_libdir "${GMP_LIBRARY}" DIRECTORY)
  set(_gmp_pc_libdirs ${PC_GMP_LIBDIR} ${PC_GMP_LIBRARY_DIRS}
      ${PC_GMP_STATIC_LIBRARY_DIRS})
  list(REMOVE_DUPLICATES _gmp_pc_libdirs)
  if(NOT _gmp_pc_libdirs)
    set(GMP_PKGCONFIG_FOUND TRUE)
  else()
    foreach(_gmp_pc_libdir IN LISTS _gmp_pc_libdirs)
      get_filename_component(_gmp_pc_libdir_abs "${_gmp_pc_libdir}" REALPATH)
      get_filename_component(_gmp_selected_libdir_abs "${_gmp_selected_libdir}" REALPATH)
      if(_gmp_pc_libdir_abs STREQUAL _gmp_selected_libdir_abs)
        set(GMP_PKGCONFIG_FOUND TRUE)
      endif()
    endforeach()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR
  VERSION_VAR GMP_VERSION)

if(GMP_FOUND)
  set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
  if(GMP_PKGCONFIG_FOUND)
    list(APPEND GMP_INCLUDE_DIRS ${PC_GMP_INCLUDE_DIRS})
  endif()
  set(GMP_LIBRARIES ${GMP_LIBRARY})
  if(NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    if(GMP_PKGCONFIG_FOUND)
      set(_gmp_interface_includes ${GMP_INCLUDE_DIR} ${PC_GMP_INCLUDE_DIRS}
          ${PC_GMP_STATIC_INCLUDE_DIRS})
      set(_gmp_interface_libs ${PC_GMP_LIBRARIES}
          ${PC_GMP_STATIC_LIBRARIES})
      set(_gmp_interface_link_dirs ${PC_GMP_LIBRARY_DIRS}
          ${PC_GMP_STATIC_LIBRARY_DIRS})
      set(_gmp_interface_compile_options ${PC_GMP_CFLAGS_OTHER}
          ${PC_GMP_STATIC_CFLAGS_OTHER})
      set(_gmp_interface_link_options ${PC_GMP_LDFLAGS_OTHER}
          ${PC_GMP_STATIC_LDFLAGS_OTHER})
    else()
      set(_gmp_interface_includes ${GMP_INCLUDE_DIR})
      set(_gmp_interface_libs)
      set(_gmp_interface_link_dirs)
      set(_gmp_interface_compile_options)
      set(_gmp_interface_link_options)
    endif()
    list(REMOVE_ITEM _gmp_interface_libs gmp)
    list(REMOVE_DUPLICATES _gmp_interface_libs)
    set_target_properties(GMP::GMP PROPERTIES
      IMPORTED_LOCATION "${GMP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${_gmp_interface_includes}"
      INTERFACE_LINK_LIBRARIES "${_gmp_interface_libs}"
      INTERFACE_LINK_DIRECTORIES "${_gmp_interface_link_dirs}"
      INTERFACE_COMPILE_OPTIONS "${_gmp_interface_compile_options}"
      INTERFACE_LINK_OPTIONS "${_gmp_interface_link_options}")
    list(APPEND GMP_INCLUDE_DIRS ${_gmp_interface_includes})
    list(APPEND GMP_LIBRARIES ${_gmp_interface_libs})
  endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
