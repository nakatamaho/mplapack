# RunStdinTest.cmake — ctest helper: run an executable with a file on stdin and
# fail if it exits non-zero or prints a failure marker.
#
# Invoked as:  cmake -DEXE=<path> -DINPUT=<file> -P RunStdinTest.cmake

if(NOT EXE OR NOT INPUT)
  message(FATAL_ERROR "RunStdinTest: EXE and INPUT must be set")
endif()

execute_process(
  COMMAND "${EXE}"
  INPUT_FILE "${INPUT}"
  RESULT_VARIABLE _rc
  OUTPUT_VARIABLE _out
  ERROR_VARIABLE _err)

message("${_out}")
if(_err)
  message("${_err}")
endif()

if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "Test driver exited with code ${_rc}")
endif()

# The LAPACK-derived drivers print human-readable pass/fail lines; treat an
# explicit failure count as a test failure.
if(_out MATCHES "[1-9][0-9]* +tests? +failed" OR _out MATCHES "\\*\\* .* failed")
  message(FATAL_ERROR "Test driver reported failures")
endif()
