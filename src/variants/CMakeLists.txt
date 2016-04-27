include_directories(../utils/htslib/htslib/
                    ../utils/bedtools/lineFileUtilities/
                    ../gtf/
                    ../utils/
                    ../utils/htslib/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/bedFile/
                    ../junctions/)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__STDC_LIMIT_MACROS")

add_library(variants
    variants_main.cc
    variants_annotator.cc)
