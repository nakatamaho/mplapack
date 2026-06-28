# Examples.  Each source is named <Name>_<backend>.cpp and uses the matching
# backend libraries.  Sources under examples/mpblas link the MPBLAS library;
# those under examples/mplapack link the MPLAPACK library (which pulls in
# MPBLAS).  The "generic" template examples and the "_opt" optimized-library
# variants from the autotools build are not part of this CMake port.

set(_mplapack_example_backends gmp mpfr qd dd double binary80 binary128)

foreach(_layer mpblas mplapack)
  foreach(_b ${_mplapack_example_backends})
    string(TOUPPER ${_b} _B)
    if(NOT MPLAPACK_ENABLE_${_B})
      continue()
    endif()

    file(GLOB_RECURSE _srcs CONFIGURE_DEPENDS
      "${CMAKE_CURRENT_SOURCE_DIR}/${_layer}/*_${_b}.cpp")

    foreach(_src ${_srcs})
      # Skip the generic/ template directory.
      if(_src MATCHES "/generic/")
        continue()
      endif()

      # Build a unique, readable target name from the path relative to examples/.
      file(RELATIVE_PATH _rel "${CMAKE_CURRENT_SOURCE_DIR}" "${_src}")
      string(REGEX REPLACE "\\.cpp$" "" _name "${_rel}")
      string(REGEX REPLACE "[/.]" "_" _name "example_${_name}")

      add_executable(${_name} "${_src}")
      target_link_libraries(${_name} PRIVATE mplapack::${_layer}_${_b})
      set_target_properties(${_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_layer}/${_b}")
    endforeach()
  endforeach()
endforeach()
