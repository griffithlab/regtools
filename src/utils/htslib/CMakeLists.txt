AUX_SOURCE_DIRECTORY(cram/ CRAM)

set(SOURCES
    bgzf.c
    faidx.c
    hts.c
    hfile.c
    hfile_net.c
    kstring.c
    knetfile.c
    sam.c)

#htslib headers have relative path to this dir
include_directories(.)

find_package(Threads)
add_library(htslib ${SOURCES} ${CRAM})
target_link_libraries (htslib ${CMAKE_THREAD_LIBS_INIT})
