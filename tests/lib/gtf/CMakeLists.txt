cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS gtf)
set(TEST_SOURCES
    test_gtf_parser.cc)

set(test_name TestGtf)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories("${PROJECT_SOURCE_DIR}/src/gtf/")
include_directories("${PROJECT_SOURCE_DIR}/src/utils/"
                    "${PROJECT_SOURCE_DIR}/src/utils/bedtools/lineFileUtilities/"
                    "${PROJECT_SOURCE_DIR}/src/utils/bedtools/gzstream/"
                    "${PROJECT_SOURCE_DIR}/src/utils/bedtools/fileType/"
                    "${PROJECT_SOURCE_DIR}/src/utils/bedtools/bedFile/")
add_executable(${test_name} ${TEST_SOURCES})
target_link_libraries(${test_name} gtest gtest_main gtf)
set(NOSTRING_FLAG "-Wno-write-strings")
set_target_properties(${test_name} PROPERTIES COMPILE_FLAGS
    ${NOSTRING_FLAG})

add_test(${TEST_SOURCES} ${test_name})
