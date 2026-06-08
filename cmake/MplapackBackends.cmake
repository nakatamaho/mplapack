# MplapackBackends.cmake
#
# Defines mplapack_add_backend(), which creates the two libraries that make up
# one precision backend:
#
#   mpblas_<backend>     (the MPBLAS sources)   -> alias mplapack::mpblas_<backend>
#   mplapack_<backend>   (the MPLAPACK sources) -> alias mplapack::mplapack_<backend>
#
# Both carry the PUBLIC compile definition that selects the backend
# (___MPLAPACK_BUILD_WITH_<BACKEND>___) so that consumers including <mplapack.h>
# automatically resolve to the matching backend header. The mplapack_<backend>
# library links the matching mpblas_<backend>.
#
# Source lists (MPBLAS_SOURCES, MPLAPACK_SOURCES) and the include directories
# are expected to be set by the top-level CMakeLists before this is used.

function(mplapack_add_backend backend macro)
  set(mpblas_tgt   mpblas_${backend})
  set(mplapack_tgt mplapack_${backend})

  add_library(${mpblas_tgt}   ${MPBLAS_SOURCES})
  add_library(${mplapack_tgt} ${MPLAPACK_SOURCES})

  add_library(mplapack::${mpblas_tgt}   ALIAS ${mpblas_tgt})
  add_library(mplapack::${mplapack_tgt} ALIAS ${mplapack_tgt})

  foreach(tgt ${mpblas_tgt} ${mplapack_tgt})
    target_compile_features(${tgt} PUBLIC cxx_std_17)
    # The backend-selecting macro must reach consumers -> PUBLIC.
    target_compile_definitions(${tgt} PUBLIC ${macro})
    target_include_directories(${tgt} PUBLIC
      "$<BUILD_INTERFACE:${MPLAPACK_INCLUDE_BUILD_DIRS}>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mplapack>")
    set_target_properties(${tgt} PROPERTIES
      VERSION ${MPLAPACK_LIB_VERSION}
      SOVERSION ${MPLAPACK_LIB_SOVERSION}
      WINDOWS_EXPORT_ALL_SYMBOLS ON)
  endforeach()

  # LAPACK layer depends on the BLAS layer of the same backend.
  target_link_libraries(${mplapack_tgt} PUBLIC ${mpblas_tgt})

  # MPLAPACK uses OpenMP in a few routines; link it when available.
  if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(${mplapack_tgt} PRIVATE OpenMP::OpenMP_CXX)
  endif()
endfunction()
