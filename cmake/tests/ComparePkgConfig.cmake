if(NOT DEFINED CMAKE_PC_DIR OR NOT DEFINED AUTOTOOLS_PC_DIR)
  message(FATAL_ERROR "CMAKE_PC_DIR and AUTOTOOLS_PC_DIR are required")
endif()

foreach(_flavor IN LISTS FLAVORS)
  set(_name "${_flavor}.pc")
  set(_cmake_pc "${CMAKE_PC_DIR}/${_name}")
  set(_autotools_pc "${AUTOTOOLS_PC_DIR}/${_name}")
  if(NOT EXISTS "${_cmake_pc}" OR NOT EXISTS "${_autotools_pc}")
    message(FATAL_ERROR
      "Missing pkg-config parity input: ${_cmake_pc} or ${_autotools_pc}")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files
            "${_cmake_pc}" "${_autotools_pc}"
    RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR
      "Autotools/CMake pkg-config mismatch for ${_name}")
  endif()
endforeach()
