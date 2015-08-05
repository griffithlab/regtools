cmake_minimum_required(VERSION 2.8)

project(regtools)


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

# Enable warnings after compiling third party dependencies.
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU|Clang")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long")
    # Uncomment when feeling persnickety.
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror")
endif()

add_subdirectory ("${PROJECT_SOURCE_DIR}/src/junctions/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/gtf/")

add_executable (regtools src/regtools.cc)
target_link_libraries (regtools junctions bedtools gtf htslib)
