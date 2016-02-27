#samtools headers have relative path to this dir
include_directories(.
                    ${CMAKE_SOURCE_DIR}/src/utils/htslib)

set(SOURCES
    bam_md.c
    bam_aux.c
    bam2bcf.c
    bam2bcf_indel.c
    bamtk.c
    bedidx.c
    errmod.c
    kprobaln.c
    sample.c
    sam_opts.c
    bam_plcmd.c)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-long-long")
add_library(samtools ${SOURCES})
