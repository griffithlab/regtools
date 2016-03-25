include_directories(../gtf/
                    ../utils/
                    ../utils/htslib/
                    ../utils/bedtools/bedFile/
                    ../utils/bedtools/lineFileUtilities/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/stringUtilities/)

add_library(junctions
    junctions_main.cc
    junctions_extractor.cc
    junctions_annotator.cc)

