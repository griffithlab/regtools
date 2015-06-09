cmake_minimum_required(VERSION 2.8)

project(regtools)
add_executable (regtools src/regtools.cc)

#versioning stuff
set (regtools_VERSION_MAJOR 0)
set (regtools_VERSION_MINOR 1)
set (regtools_VERSION_PATCH 0)
configure_file (
    "${PROJECT_SOURCE_DIR}/src/version.h.in"
    "${PROJECT_BINARY_DIR}/version.h"
)
include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/bedtools/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/samtools/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/junctions/")
target_link_libraries (regtools junctions bedtools samtools)

