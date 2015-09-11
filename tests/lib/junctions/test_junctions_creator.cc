#include <gtest/gtest.h>
#include <sstream>
#include "junctions_creator.h"

class JunctionsCreateTest : public ::testing::Test {
    public:
        JunctionsCreator j1;
};

TEST_F(JunctionsCreateTest, ParseInput) {
    int argc = 2;
    char * argv[] = {"regtools", "test_input.bam"};
    int ret = j1.parse_options(argc, argv);
    string expected_bam("test_input.bam");
    ASSERT_EQ(expected_bam, j1.get_bam());
    ASSERT_EQ(0, ret);
}

TEST_F(JunctionsCreateTest, ParseEmptyInput) {
    int argc = 1;
    char * argv[] = {"regtools"};
    int ret = j1.parse_options(argc, argv);
    ASSERT_EQ(1, ret) << "Error parsing as expected";
}

TEST_F(JunctionsCreateTest, Usage) {
    ostringstream out, out2;
    out << "\nUsage:\t\t" << "regtools junctions create [options] indexed_alignments.bam";
    out << "\nOptions:";
    out << "\t" << "-a FILE\tMinimum anchor length. Junctions which satisfy a minimum "
                     "anchor length on both sides are reported. [8]";
    out << "\n\t\t" << "-i FILE\tMinimum intron length. [70]";
    out << "\n\t\t" << "-I FILE\tMaximum intron length. [500000]";
    out << "\n\t\t" << "-o FILE\tThe file to write output to. [STDOUT]";
    out << "\n\t\t" << "-r STR\tThe region to identify junctions "
                     "in \"chr:start-end\" format. Entire BAM by default.";
    out << "\n";
    j1.usage(out2);
    ASSERT_EQ(out.str(), out2.str()) << "Error parsing as expected";
}
