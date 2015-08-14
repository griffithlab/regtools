cmake_minimum_required(VERSION 2.8)

#Add build-common cmake modules path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build-common/cmake)
#include(ProjectHelper)
include(TestHelper)

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
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/htslib/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/gtest-1.7.0/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/junctions/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/gtf/")
target_link_libraries (regtools junctions bedtools gtf htslib)

#Testing
enable_testing()
add_subdirectory(tests/lib) #unit-tests
add_subdirectory(tests/integration-test) #integration-tests
