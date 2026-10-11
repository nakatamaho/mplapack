# Copy installed .pc files, pointing their location variables at the
# directory the test suite installed into.
#
# The installed files carry CMAKE_INSTALL_PREFIX from configure time (they
# must stay identical to the Autotools ones), while the consume tests install
# with "cmake --install --prefix".  Only the prefix, exec_prefix, libdir and
# includedir lines are rewritten, and only where they start with FROM; every
# other line, including dependency flags, is copied unchanged.
#
#   cmake -DSRC=<dir> -DDST=<dir> -DFROM=<prefix> -DTO=<prefix> -P RelocatePkgConfig.cmake

foreach(_var SRC DST FROM TO)
  if(NOT DEFINED ${_var})
    message(FATAL_ERROR "${_var} is required")
  endif()
endforeach()

file(REMOVE_RECURSE "${DST}")
file(MAKE_DIRECTORY "${DST}")
file(GLOB _pc_files "${SRC}/*.pc")
if(NOT _pc_files)
  message(FATAL_ERROR "No .pc files in ${SRC}")
endif()

string(LENGTH "${FROM}" _from_length)
foreach(_pc IN LISTS _pc_files)
  get_filename_component(_name "${_pc}" NAME)
  file(STRINGS "${_pc}" _lines)
  set(_out "")
  foreach(_line IN LISTS _lines)
    if(_line MATCHES "^(prefix|exec_prefix|libdir|includedir)=(.*)$")
      set(_key "${CMAKE_MATCH_1}")
      set(_value "${CMAKE_MATCH_2}")
      string(FIND "${_value}" "${FROM}" _pos)
      if(_pos EQUAL 0)
        string(SUBSTRING "${_value}" ${_from_length} -1 _rest)
        set(_line "${_key}=${TO}${_rest}")
      endif()
    endif()
    string(APPEND _out "${_line}\n")
  endforeach()
  file(WRITE "${DST}/${_name}" "${_out}")
endforeach()
