include_directories(../gtf/
                    ../utils/
                    ../utils/bedtools/bedFile/
                    ../utils/bedtools/lineFileUtilities/
                    ../utils/bedtools/gzstream/
                    ../utils/bedtools/fileType/
                    ../utils/bedtools/stringUtilities/
                    ../utils/htslib/
                    ../utils/samtools/
                    ../variants/
                    ../junctions/)

add_library(cis-ase
    cis_ase_identifier.cc
    cis_ase_main.cc
    bam_plcmd_regtools.c)
