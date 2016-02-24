AUX_SOURCE_DIRECTORY(cram/ CRAM)

set(SOURCES
    bgzf.c
    faidx.c
    hts.c
    hfile.c
    hfile_net.c
    kfunc.c
    kstring.c
    knetfile.c
    md5.c
    tbx.c
    vcf.c
    sam.c)

#htslib headers have relative path to this dir
include_directories(.)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-long-long")
find_package(Threads)
add_library(htslib ${SOURCES} ${CRAM})
target_link_libraries (htslib ${CMAKE_THREAD_LIBS_INIT})
