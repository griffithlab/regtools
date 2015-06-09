include_directories(bedFile/
                    Fasta/
                    fileType/
                    gzstream/
                    lineFileUtilities/
                    stringUtilities/)
set(SOURCES
    gzstream/gzstream.C
    fileType/fileType.cpp
    bedFile/bedFile.cpp)

add_library(bedtools ${SOURCES})

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
        include_directories(${ZLIB_INCLUDE_DIRS})
        target_link_libraries(bedtools ${ZLIB_LIBRARIES} )
endif(ZLIB_FOUND)

