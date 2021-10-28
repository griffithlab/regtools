cmake_minimum_required(VERSION 2.8)

project(regtools)

#Add build-common cmake modules path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build-common/cmake)
include(ProjectHelper)
include(TestHelper)



#versioning stuff
set (regtools_VERSION_MAJOR 0)
set (regtools_VERSION_MINOR 0)
set (regtools_VERSION_PATCH 1)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/version.h.in"
    "${PROJECT_BINARY_DIR}/version.h"
)
include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/bedtools/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/htslib/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/samtools/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/gtest-1.7.0/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/utils/rmath/")

# Enable warnings after compiling third party dependencies.
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU|Clang")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-long-long")
    # Uncomment when feeling persnickety.
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror")
endif()

add_subdirectory ("${PROJECT_SOURCE_DIR}/src/gtf/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/junctions/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/cis-ase/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/cis-splice-effects/")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/variants/")

#The main executable
add_executable (regtools src/regtools.cc)
target_link_libraries (regtools junctions variants
                       cis-ase bedtools gtf rmath samtools htslib cis-splice-effects )

#Testing
enable_testing()
add_subdirectory(tests/lib) #unit-tests
add_subdirectory(tests/integration-test) #integration-tests
