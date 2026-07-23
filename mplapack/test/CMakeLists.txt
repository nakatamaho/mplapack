# Test suite (LAPACK-derived lin/ and eig/ drivers), ported from the autotools
# test harness.  For each enabled backend this builds:
#
#   libmatgen_<b>        — test matrix generators (mplapack/test/matgen)
#   liblin_<b>           — linear-equation test routines (lin/common)
#   liblin_override_<b>  — test versions of Mxerbla/iMlaenv/Mxlaenv
#   xlintst{R,C,rfR,rfC}_<b>  — the lin test drivers
#   libeig_<b> / libeig_override_<b> / xeigtst{R,C}_<b> / xdmdeigtst{R,C}_<b>
#
# and registers ctest cases that feed each driver its .in dataset on stdin.
#
# The override library provides alternative definitions of a few symbols that
# also live in libmplapack; it is linked first and -Wl,--allow-multiple-definition
# makes those definitions win (mirroring the autotools whole-archive trick).

set(_test_dir       "${CMAKE_CURRENT_SOURCE_DIR}")
set(_matgen_dir     "${_test_dir}/matgen")
set(_lin_common     "${_test_dir}/lin/common")
set(_eig_common     "${_test_dir}/eig/common")
set(_fable_dir      "${CMAKE_SOURCE_DIR}/fable")

set(_override_files Mxerbla iMlaenv Mxlaenv)
set(_lin_mains      Rchkaa Cchkaa Rchkrfp Cchkrfp)
set(_eig_mains      Rchkee Cchkee Rchkdmd Cchkdmd)

# MPLAPACK_LD_ALLOW_MULTIDEF is probed once by the top-level CMakeLists.

# matgen / lin-library / eig-library source lists (shared across backends).
file(GLOB _matgen_srcs CONFIGURE_DEPENDS "${_matgen_dir}/*.cpp")

file(GLOB _lin_lib_srcs CONFIGURE_DEPENDS "${_lin_common}/*.cpp")
foreach(_x ${_override_files} ${_lin_mains})
  list(FILTER _lin_lib_srcs EXCLUDE REGEX "/${_x}\\.cpp$")
endforeach()

file(GLOB _eig_lib_srcs CONFIGURE_DEPENDS "${_eig_common}/*.cpp")
foreach(_x ${_override_files} ${_eig_mains})
  list(FILTER _eig_lib_srcs EXCLUDE REGEX "/${_x}\\.cpp$")
endforeach()

# Helper: apply the common test compile settings to a target.
function(_mplapack_test_setup tgt backend_target internal)
  # PRIVATE: these are static archives; they only need the backend's headers,
  # the ___MPLAPACK_BUILD_WITH_<backend>___ macro and the fable headers to
  # compile.  We deliberately do NOT propagate the backend link to consumers --
  # the driver executables link the backend explicitly inside a LINK_GROUP so
  # that ordering/override resolution is under our control.
  target_link_libraries(${tgt} PRIVATE ${backend_target})
  target_include_directories(${tgt} PRIVATE "${_fable_dir}")
  if(internal)
    target_compile_definitions(${tgt} PRIVATE ___MPLAPACK_INTERNAL___)
  endif()
  if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(${tgt} PRIVATE OpenMP::OpenMP_CXX)
  endif()
endfunction()

# lin .in datasets -> driver suffix.
set(_lin_R_inputs    Rtest)        # -> xlintstR
set(_lin_C_inputs    Ctest)        # -> xlintstC
set(_lin_rfR_inputs  Rtest_rfp)    # -> xlintstrfR
set(_lin_rfC_inputs  Ctest_rfp)    # -> xlintstrfC

# eig .in datasets per driver (from the autotools check-local recipes).
set(_eig_R_inputs   nep sep se2 svd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb glm gqr gsv csd lse)
set(_eig_C_inputs   nep sep se2 svd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb glm gqr gsv csd lse)
set(_dmd_R_inputs   Rdmd)
set(_dmd_C_inputs   Cdmd)

function(_mplapack_register_stdin_test name exe input_file)
  if(NOT EXISTS "${input_file}")
    return()
  endif()
  add_test(NAME ${name}
    COMMAND ${CMAKE_COMMAND}
      -DEXE=$<TARGET_FILE:${exe}>
      -DINPUT=${input_file}
      -P "${CMAKE_SOURCE_DIR}/cmake/RunStdinTest.cmake")
endfunction()

set(_backends gmp mpfr qd dd double binary80 binary128)
foreach(_b ${_backends})
  string(TOUPPER ${_b} _B)
  if(NOT MPLAPACK_ENABLE_${_B})
    continue()
  endif()

  set(_backend_libraries mplapack_${_b})
  if(TARGET mplapack_${_b}_opt)
    list(APPEND _backend_libraries mplapack_${_b}_opt)
  endif()

  foreach(_backend_library IN LISTS _backend_libraries)
    if(_backend_library MATCHES "_opt$")
      set(_variant "${_b}_opt")
      set(_register_stdin_tests FALSE)
    else()
      set(_variant "${_b}")
      set(_register_stdin_tests TRUE)
    endif()

    # --- libraries --------------------------------------------------------
    add_library(matgen_${_variant}        STATIC ${_matgen_srcs})
    add_library(lin_${_variant}           STATIC ${_lin_lib_srcs})
    add_library(lin_override_${_variant}  STATIC)
    add_library(eig_${_variant}           STATIC ${_eig_lib_srcs})
    add_library(eig_override_${_variant}  STATIC)
    foreach(_o ${_override_files})
      target_sources(lin_override_${_variant} PRIVATE "${_lin_common}/${_o}.cpp")
      target_sources(eig_override_${_variant} PRIVATE "${_eig_common}/${_o}.cpp")
    endforeach()

    _mplapack_test_setup(matgen_${_variant}       ${_backend_library} TRUE)
    _mplapack_test_setup(lin_${_variant}          ${_backend_library} TRUE)
    _mplapack_test_setup(lin_override_${_variant} ${_backend_library} FALSE)
    _mplapack_test_setup(eig_${_variant}          ${_backend_library} TRUE)
    _mplapack_test_setup(eig_override_${_variant} ${_backend_library} FALSE)

    # --- drivers ----------------------------------------------------------
    # name : driver-suffix : main-source : dataset-list-variable
    set(_lin_drivers
      "xlintstR_${_variant}|${_lin_common}/Rchkaa.cpp|_lin_R_inputs"
      "xlintstC_${_variant}|${_lin_common}/Cchkaa.cpp|_lin_C_inputs"
      "xlintstrfR_${_variant}|${_lin_common}/Rchkrfp.cpp|_lin_rfR_inputs"
      "xlintstrfC_${_variant}|${_lin_common}/Cchkrfp.cpp|_lin_rfC_inputs")
    set(_eig_drivers
      "xeigtstR_${_variant}|${_eig_common}/Rchkee.cpp|_eig_R_inputs"
      "xeigtstC_${_variant}|${_eig_common}/Cchkee.cpp|_eig_C_inputs"
      "xdmdeigtstR_${_variant}|${_eig_common}/Rchkdmd.cpp|_dmd_R_inputs"
      "xdmdeigtstC_${_variant}|${_eig_common}/Cchkdmd.cpp|_dmd_C_inputs")

    foreach(_spec ${_lin_drivers} ${_eig_drivers})
      string(REPLACE "|" ";" _parts "${_spec}")
      list(GET _parts 0 _exe)
      list(GET _parts 1 _main)
      list(GET _parts 2 _inputs_var)

      add_executable(${_exe} "${_main}")
      target_compile_definitions(${_exe} PRIVATE ___MPLAPACK_INTERNAL___)
      target_include_directories(${_exe} PRIVATE "${_fable_dir}")
      if(_b STREQUAL "dd" AND MPLAPACK_HAS_FP_CONTRACT_OFF)
        target_compile_options(${_exe} PRIVATE -ffp-contract=off)
      endif()
      if(_exe MATCHES "^xeig|^xdmd")
        set(_kind eig)
      else()
        set(_kind lin)
      endif()
      # The override lib provides the test versions of Mxerbla, iMlaenv and
      # Mxlaenv (and defines srnamt/infot) that must REPLACE the real ones in
      # the backend libraries -- e.g. the error-exit tests rely on the
      # override Mxerbla recording the error instead of printing and aborting.
      # Linking the override archive with --whole-archive forces all its
      # objects in first, so its definitions win over the duplicates in the
      # merged backend library (kept first by --allow-multiple-definition).
      # The remaining static archives have circular dependencies, so wrap
      # them in a RESCAN group (--start-group/--end-group).
      target_link_libraries(${_exe} PRIVATE
        "$<LINK_LIBRARY:WHOLE_ARCHIVE,${_kind}_override_${_variant}>"
        "$<LINK_GROUP:RESCAN,${_kind}_${_variant},matgen_${_variant},${_backend_library}>")
      if(MPLAPACK_LD_ALLOW_MULTIDEF)
        target_link_options(${_exe} PRIVATE "-Wl,--allow-multiple-definition")
      endif()
      if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(${_exe} PRIVATE OpenMP::OpenMP_CXX)
      endif()

      if(_kind STREQUAL lin)
        set(_dataset_dir "${_test_dir}/lin")
      else()
        set(_dataset_dir "${_test_dir}/eig")
      endif()
      if(_register_stdin_tests)
        foreach(_in ${${_inputs_var}})
          # ${_exe} embeds the backend flavor and R/C/rf/dmd driver variant, so
          # it disambiguates datasets shared between drivers (e.g. nep, svd).
          set(_input_file "${_dataset_dir}/${_in}.in")
          _mplapack_register_stdin_test(
            "${_exe}_${_in}" ${_exe} "${_input_file}")
        endforeach()
      endif()
    endforeach()
  endforeach()
endforeach()
