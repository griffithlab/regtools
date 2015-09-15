#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include "junctions_creator.h"

class JunctionsCreateTest : public ::testing::Test {
    public:
        JunctionsCreator j1;
};

TEST_F(JunctionsCreateTest, ParseInput) {
    int argc = 2;
    char * argv[] = {"create", "test_input.bam"};
    int ret = j1.parse_options(argc, argv);
    string expected_bam("test_input.bam");
    ASSERT_EQ(expected_bam, j1.get_bam());
    ASSERT_EQ(0, ret);
}

TEST_F(JunctionsCreateTest, ParseNoInput) {
    int argc = 1;
    char * argv[] = {"create"};
    ASSERT_THROW(j1.parse_options(argc, argv), std::runtime_error);
}

TEST_F(JunctionsCreateTest, ParseIncorrectOption) {
    int argc = 2;
    char * argv[] = {"create", "-k", "24", "test_input.bam"};
    ASSERT_THROW(j1.parse_options(argc, argv), std::runtime_error);
}

TEST_F(JunctionsCreateTest, Usage) {
    ostringstream out, out2;
    out << "\nUsage:\t\t" << "regtools junctions create [options] indexed_alignments.bam";
    out << "\nOptions:";
    out << "\t" << "-a INT\tMinimum anchor length. Junctions which satisfy a minimum "
                     "anchor length on both sides are reported. [8]";
    out << "\n\t\t" << "-i INT\tMinimum intron length. [70]";
    out << "\n\t\t" << "-I INT\tMaximum intron length. [500000]";
    out << "\n\t\t" << "-o FILE\tThe file to write output to. [STDOUT]";
    out << "\n\t\t" << "-r STR\tThe region to identify junctions "
                     "in \"chr:start-end\" format. Entire BAM by default.";
    out << "\n";
    j1.usage(out2);
    ASSERT_EQ(out.str(), out2.str()) << "Error parsing as expected";
}
