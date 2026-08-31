cmake_minimum_required(VERSION 3.20)
project(mplapack_consume LANGUAGES CXX)

if(NOT MPLAPACK_CONSUME_TARGET)
  message(FATAL_ERROR "MPLAPACK_CONSUME_TARGET is required")
endif()
if(NOT MPLAPACK_CONSUME_SOURCE)
  message(FATAL_ERROR "MPLAPACK_CONSUME_SOURCE is required")
endif()

find_package(mplapack CONFIG REQUIRED)

if(NOT TARGET "mplapack::${MPLAPACK_CONSUME_TARGET}")
  message(FATAL_ERROR
    "Package does not export mplapack::${MPLAPACK_CONSUME_TARGET}")
endif()

add_executable(mplapack_consume "${MPLAPACK_CONSUME_SOURCE}")
target_link_libraries(mplapack_consume PRIVATE
  "mplapack::${MPLAPACK_CONSUME_TARGET}")
