# Apply the bundled LAPACK patches to an extracted reference LAPACK tree, in
# the order external/lapack/Makefile.am applies them.
#
#   cmake -DSOURCE_DIR=<lapack tree> -DPATCH_DIR=<patches/<version>>
#         -DPATCHES=<a|b|...> -DPATCH=<patch program> -P ApplyLapackPatches.cmake

foreach(_var SOURCE_DIR PATCH_DIR PATCHES PATCH)
  if(NOT DEFINED ${_var})
    message(FATAL_ERROR "${_var} is required")
  endif()
endforeach()

string(REPLACE "|" ";" _patches "${PATCHES}")
foreach(_patch IN LISTS _patches)
  execute_process(
    COMMAND "${PATCH}" -p1 -N -i "${PATCH_DIR}/${_patch}"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _output)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Applying ${_patch} failed:\n${_output}")
  endif()
endforeach()
