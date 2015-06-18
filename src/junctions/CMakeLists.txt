include_directories(../gtf/
                    ../utils/
                    ../utils/bedtools/bedFile/
                    ../utils/bedtools/lineFileUtilities/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/stringUtilities/
                    ../utils/samtools/)

add_library(junctions
    junctions_main.cc
    junctions_annotator.cc)

