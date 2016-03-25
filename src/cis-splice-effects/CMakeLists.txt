include_directories(../gtf/
                    ../utils/
                    ../utils/bedtools/bedFile/
                    ../utils/bedtools/lineFileUtilities/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/stringUtilities/
                    ../utils/htslib/
                    ../variants/
                    ../junctions/)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__STDC_LIMIT_MACROS")

add_library(cis-splice-effects
    cis_splice_effects_identifier.cc
    cis_splice_effects_main.cc)
