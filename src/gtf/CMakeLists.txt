include_directories(../utils/
                    ../utils/bedtools/lineFileUtilities/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/bedFile/)
add_library(gtf
    gtf_parser.cc)

