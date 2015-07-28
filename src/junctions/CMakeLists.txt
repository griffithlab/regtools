include_directories(../gtf/
                    ../utils/
                    ../utils/bedtools/bedFile/
                    ../utils/bedtools/lineFileUtilities/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/stringUtilities/
                    ../utils/htslib/htslib/)

add_library(junctions
    junctions_main.cc
    junctions_creator.cc
    junctions_annotator.cc)

