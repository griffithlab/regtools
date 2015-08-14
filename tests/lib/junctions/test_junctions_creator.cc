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
    ostringstream out1, out2;
    out1 << "\nUsage:\t\t" << "regtools junctions create [options] alignments.bam";
    out1 << "\n";
    j1.usage(out2);
    ASSERT_EQ(out1.str(), out2.str()) << "Error parsing as expected";
}
