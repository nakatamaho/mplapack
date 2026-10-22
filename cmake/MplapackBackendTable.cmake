# MplapackBackendTable.cmake — read backends.txt, the MPLAPACK backend table.
#
# Defines
#   MPLAPACK_BACKENDS                backend names in table order
#   MPLAPACK_BACKEND_<name>_REAL     C++ type of REAL
#   MPLAPACK_BACKEND_<name>_COMPLEX  C++ type of COMPLEX
#   MPLAPACK_BACKEND_<name>_DEFAULT  default of MPLAPACK_ENABLE_<NAME>
#   MPLAPACK_BACKEND_<name>_TRAITS   list of build traits
# and the functions below.  backends.txt documents the columns and traits.

set(MPLAPACK_BACKEND_TABLE "${CMAKE_CURRENT_LIST_DIR}/../backends.txt")
get_filename_component(MPLAPACK_BACKEND_TABLE "${MPLAPACK_BACKEND_TABLE}" ABSOLUTE)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${MPLAPACK_BACKEND_TABLE}")

# Traits handled by mplapack_configure_backend_traits (MplapackBackends.cmake).
set(_mplapack_known_backend_traits gmp mpfr qd nofma binary128libs)

set(MPLAPACK_BACKENDS "")
file(STRINGS "${MPLAPACK_BACKEND_TABLE}" _mplapack_backend_rows
  REGEX "^[^# \t]")
foreach(_row IN LISTS _mplapack_backend_rows)
  string(REGEX REPLACE "[ \t]+" ";" _fields "${_row}")
  list(LENGTH _fields _nfields)
  if(NOT _nfields EQUAL 5)
    message(FATAL_ERROR "${MPLAPACK_BACKEND_TABLE}: expected 5 columns: ${_row}")
  endif()
  list(GET _fields 0 _name)
  if(NOT _name MATCHES "^[a-z][a-z0-9]*$")
    message(FATAL_ERROR "${MPLAPACK_BACKEND_TABLE}: bad backend name '${_name}'")
  endif()
  if(_name IN_LIST MPLAPACK_BACKENDS)
    message(FATAL_ERROR "${MPLAPACK_BACKEND_TABLE}: duplicate backend '${_name}'")
  endif()
  list(APPEND MPLAPACK_BACKENDS "${_name}")
  list(GET _fields 1 MPLAPACK_BACKEND_${_name}_REAL)
  list(GET _fields 2 MPLAPACK_BACKEND_${_name}_COMPLEX)
  list(GET _fields 3 MPLAPACK_BACKEND_${_name}_DEFAULT)
  list(GET _fields 4 _traits)
  if(_traits STREQUAL "-")
    set(MPLAPACK_BACKEND_${_name}_TRAITS "")
  else()
    string(REPLACE "," ";" MPLAPACK_BACKEND_${_name}_TRAITS "${_traits}")
  endif()
  foreach(_trait IN LISTS MPLAPACK_BACKEND_${_name}_TRAITS)
    if(NOT _trait IN_LIST _mplapack_known_backend_traits)
      message(FATAL_ERROR
        "${MPLAPACK_BACKEND_TABLE}: unknown trait '${_trait}' for ${_name}")
    endif()
  endforeach()
endforeach()

# Traits of every backend whose MPLAPACK_ENABLE_<NAME> option is on.
function(mplapack_enabled_backend_traits out)
  set(_all "")
  foreach(_b IN LISTS MPLAPACK_BACKENDS)
    string(TOUPPER "${_b}" _B)
    if(MPLAPACK_ENABLE_${_B})
      list(APPEND _all ${MPLAPACK_BACKEND_${_b}_TRAITS})
    endif()
  endforeach()
  list(REMOVE_DUPLICATES _all)
  set(${out} "${_all}" PARENT_SCOPE)
endfunction()

# Backend that a library target belongs to: mplapack_<name>[_<variant>].
function(mplapack_target_backend out target)
  string(REGEX REPLACE "^mplapack_([a-z0-9]+).*$" "\\1" _b "${target}")
  if(NOT _b IN_LIST MPLAPACK_BACKENDS)
    set(_b "")
  endif()
  set(${out} "${_b}" PARENT_SCOPE)
endfunction()
