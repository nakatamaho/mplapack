# CMake package-consumption checks for build-tree and installed exports.

add_test(NAME source_manifests
  COMMAND "${PROJECT_SOURCE_DIR}/misc/check_source_manifests.sh"
          "${PROJECT_SOURCE_DIR}")

set(_consume_source_dir "${CMAKE_CURRENT_SOURCE_DIR}/consume")
set(_consume_root "${CMAKE_CURRENT_BINARY_DIR}/consume")
set(_consume_install_prefix "${_consume_root}/install")
set(_consume_install_fixture mplapack_consume_install)

# A child configure does not inherit dependency cache entries from this build.
# Forward the locations already resolved by the parent so a non-system QD
# installation works for both the build-tree and installed-package checks.
set(_consume_dependency_args "")
foreach(_var QD_INCLUDE_DIR QD_LIBRARY
             GMP_INCLUDE_DIR GMP_LIBRARY
             MPFR_INCLUDE_DIR MPFR_LIBRARY
             MPC_INCLUDE_DIR MPC_LIBRARY)
  if(DEFINED ${_var} AND NOT "${${_var}}" STREQUAL "")
    list(APPEND _consume_dependency_args "-D${_var}=${${_var}}")
  endif()
endforeach()

add_test(NAME consume_install_tree
  COMMAND "${CMAKE_COMMAND}" --install "${CMAKE_BINARY_DIR}"
          --prefix "${_consume_install_prefix}" --config $<CONFIG>)
set_tests_properties(consume_install_tree PROPERTIES
  FIXTURES_SETUP "${_consume_install_fixture}"
  RESOURCE_LOCK mplapack_consume_install)

function(mplapack_add_consume_tests flavor)
  set(_target "mplapack_${flavor}")
  if(NOT TARGET "${_target}")
    return()
  endif()

  set(_buildtree_dir "${_consume_root}/buildtree_${flavor}")
  set(_buildtree_fixture "mplapack_consume_buildtree_${flavor}")
  if(flavor MATCHES "^mpfr")
    set(_consume_source "${CMAKE_CURRENT_SOURCE_DIR}/consume/mpfr_main.cpp")
  else()
    set(_consume_source "${CMAKE_CURRENT_SOURCE_DIR}/consume/main.cpp")
  endif()
  add_test(NAME "consume_buildtree_${flavor}_configure"
    COMMAND "${CMAKE_COMMAND}"
      -S "${_consume_source_dir}"
      -B "${_buildtree_dir}"
      "-Dmplapack_DIR=${CMAKE_BINARY_DIR}"
      "-DMPLAPACK_CONSUME_TARGET=${_target}"
      "-DMPLAPACK_CONSUME_SOURCE=${_consume_source}"
      ${_consume_dependency_args})
  set_tests_properties("consume_buildtree_${flavor}_configure" PROPERTIES
    FIXTURES_SETUP "${_buildtree_fixture}")
  add_test(NAME "consume_buildtree_${flavor}_build"
    COMMAND "${CMAKE_COMMAND}" --build "${_buildtree_dir}" --config $<CONFIG>)
  set_tests_properties("consume_buildtree_${flavor}_build" PROPERTIES
    FIXTURES_REQUIRED "${_buildtree_fixture}")

  set(_installed_dir "${_consume_root}/installed_${flavor}")
  set(_installed_fixture "mplapack_consume_installed_${flavor}")
  add_test(NAME "consume_installed_${flavor}_configure"
    COMMAND "${CMAKE_COMMAND}"
      -S "${_consume_source_dir}"
      -B "${_installed_dir}"
      "-Dmplapack_DIR=${_consume_install_prefix}/${CMAKE_INSTALL_LIBDIR}/cmake/mplapack"
      "-DMPLAPACK_CONSUME_TARGET=${_target}"
      "-DMPLAPACK_CONSUME_SOURCE=${_consume_source}"
      ${_consume_dependency_args})
  set_tests_properties("consume_installed_${flavor}_configure" PROPERTIES
    FIXTURES_REQUIRED "${_consume_install_fixture}"
    FIXTURES_SETUP "${_installed_fixture}"
    RESOURCE_LOCK mplapack_consume_install)
  add_test(NAME "consume_installed_${flavor}_build"
    COMMAND "${CMAKE_COMMAND}" --build "${_installed_dir}" --config $<CONFIG>)
  set_tests_properties("consume_installed_${flavor}_build" PROPERTIES
    FIXTURES_REQUIRED "${_installed_fixture}")
endfunction()

mplapack_add_consume_tests(dd)
mplapack_add_consume_tests(dd_opt)
mplapack_add_consume_tests(mpfr)
mplapack_add_consume_tests(mpfr_opt)

# Exercise the original underlinking failure in a fresh process when this is
# a shared CMake build.  The dependency directories are added only to the
# test environment; they are not embedded as an RPATH in the installed target.
if(BUILD_SHARED_LIBS AND TARGET mplapack_mpfr)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    set(_mpfr_runtime_dirs "$<TARGET_FILE_DIR:mplapack_mpfr>")
    foreach(_dependency_library MPC_LIBRARY MPFR_LIBRARY GMP_LIBRARY)
      if(DEFINED ${_dependency_library} AND
         NOT "${${_dependency_library}}" STREQUAL "")
        get_filename_component(_dependency_dir
          "${${_dependency_library}}" DIRECTORY)
        list(APPEND _mpfr_runtime_dirs "${_dependency_dir}")
      endif()
    endforeach()
    string(JOIN ":" _mpfr_runtime_path ${_mpfr_runtime_dirs})
    add_test(NAME mpfr_shared_clean_load
      COMMAND "${CMAKE_COMMAND}" -E env
        "LD_LIBRARY_PATH=${_mpfr_runtime_path}"
        "DYLD_LIBRARY_PATH=${_mpfr_runtime_path}"
        "PATH=${_mpfr_runtime_path}:$ENV{PATH}"
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/misc/check_mpfr_shared_load.py"
        "$<TARGET_FILE:mplapack_mpfr>")
  endif()
  if(UNIX AND NOT APPLE)
    foreach(_mpfr_target mplapack_mpfr mplapack_mpfr_opt)
      if(TARGET ${_mpfr_target})
        add_test(NAME "${_mpfr_target}_elf_dependencies"
          COMMAND "${PROJECT_SOURCE_DIR}/misc/check_mpfr_elf_dependencies.sh"
                  "$<TARGET_FILE:${_mpfr_target}>")
      endif()
    endforeach()
  endif()
endif()

set(MPLAPACK_AUTOTOOLS_PKGCONFIG_DIR "" CACHE PATH
    "Autotools build directory containing per-flavor .pc files to compare")
if(MPLAPACK_AUTOTOOLS_PKGCONFIG_DIR)
  add_test(NAME pkgconfig_autotools_cmake_parity
    COMMAND "${CMAKE_COMMAND}"
      "-DCMAKE_PC_DIR=${CMAKE_BINARY_DIR}/pkgconfig"
      "-DAUTOTOOLS_PC_DIR=${MPLAPACK_AUTOTOOLS_PKGCONFIG_DIR}"
      "-DFLAVORS=${MPLAPACK_INSTALL_TARGETS}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/ComparePkgConfig.cmake")
endif()
