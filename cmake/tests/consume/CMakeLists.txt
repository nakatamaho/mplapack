cmake_minimum_required(VERSION 3.20)
project(mplapack_consume LANGUAGES CXX)

if(NOT MPLAPACK_CONSUME_TARGET)
  message(FATAL_ERROR "MPLAPACK_CONSUME_TARGET is required")
endif()

find_package(mplapack CONFIG REQUIRED)

if(NOT TARGET "mplapack::${MPLAPACK_CONSUME_TARGET}")
  message(FATAL_ERROR
    "Package does not export mplapack::${MPLAPACK_CONSUME_TARGET}")
endif()

if(NOT MPLAPACK_CONSUME_SOURCE)
  set(MPLAPACK_CONSUME_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
endif()
add_executable(mplapack_consume "${MPLAPACK_CONSUME_SOURCE}")
target_link_libraries(mplapack_consume PRIVATE
  "mplapack::${MPLAPACK_CONSUME_TARGET}")
