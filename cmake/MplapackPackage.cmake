# MplapackPackage.cmake — CPack configuration.
#
# Provides:
#   cpack --config CPackSourceConfig.cmake   (or: cmake --build build -t package_source)
#     -> a source tarball, the rough equivalent of `make dist`.
#   cpack                                     (or: cmake --build build -t package)
#     -> a binary archive of the install tree.

set(CPACK_PACKAGE_NAME "mplapack")
set(CPACK_PACKAGE_VENDOR "Nakata Maho")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "${PROJECT_HOMEPAGE_URL}")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_VERBATIM_VARIABLES ON)

# Binary archive of the installed tree.
set(CPACK_GENERATOR "TXZ")

# Source tarball.
set(CPACK_SOURCE_GENERATOR "TXZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "mplapack-${PROJECT_VERSION}")

# Keep VCS, build trees and generated build-system output out of the source
# package (mirrors the .gitignore set, since CPack source packs the source dir).
set(CPACK_SOURCE_IGNORE_FILES
  "/\\.git/"
  "/\\.github/"
  "/build.*/"
  "${CMAKE_BINARY_DIR}/"
  # autotools-generated output
  "/autom4te\\.cache/"
  "/configure$"
  "/aclocal\\.m4$"
  "/config\\.(log|status|guess|sub)$"
  "/(compile|depcomp|install-sh|missing|libtool|ltmain\\.sh|test-driver|INSTALL)$"
  "/Makefile$"
  "/Makefile\\.in$"
  "/m4/(libtool|ltoptions|ltsugar|ltversion|lt~obsolete)\\.m4$"
  "/include/mplapack_config\\.h(\\.in)?$"
  "/mplapack\\.pc$"
  # build artifacts
  "\\.(o|lo|la|a|so)$"
  "/\\.libs/"
  "/\\.deps/"
  "\\.dirstamp$")

include(CPack)
