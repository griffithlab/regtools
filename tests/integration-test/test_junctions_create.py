#!/usr/bin/env python

from integrationtest import IntegrationTest, main
import unittest

class TestCreate(IntegrationTest, unittest.TestCase):
    def test_junctions_create_anchor(self):
        bam1 = self.inputFiles("test_hcc1395.bam")[0]
        output_file = self.tempFile("create.out")
        print "BAM1 is ", bam1
        for anchor in ["", "30"]:
            expected_file = self.inputFiles("junctions-create/expected-a" +
                                            anchor + ".out")[0]
            if anchor != "":
                anchor = "-a " + anchor
            params = ["junctions", "create", anchor, "-o", output_file, bam1]
            rv, err = self.execute(params)
            self.assertEqual(rv, 0)
            self.assertFilesEqual(expected_file, output_file)

    def test_junctions_create_intron_size(self):
        bam1 = self.inputFiles("test_hcc1395.bam")[0]
        output_file = self.tempFile("create.out")
        min_intron = "8039"
        max_intron = "8039"
        expected_file = self.inputFiles("junctions-create/expected-i" +
                min_intron + "-I" + max_intron +
                ".out")[0]
        params = ["junctions", "create", "-o", output_file,
                  "-i", min_intron, "-I", max_intron, bam1]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)
        self.assertFilesEqual(expected_file, output_file)

    def test_region(self):
        bam1 = self.inputFiles("test_hcc1395.bam")[0]
        output_file = self.tempFile("create.out")
        region = "1:22405013-22405020"
        expected_file = self.inputFiles("junctions-create/expected-r" +
                region + ".out")[0]
        params = ["junctions", "create", "-o", output_file, "-r", region,
                  bam1]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)
        self.assertFilesEqual(expected_file, output_file)

    def test_no_bam(self):
        output_file = self.tempFile("create.out")
        params = ["junctions", "create", "-o", output_file]
        rv, err = self.execute(params)
        self.assertEqual(rv, 1)

    def test_help(self):
        output_file = self.tempFile("create.out")
        params = ["junctions", "create", "-h"]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)

if __name__ == "__main__":
    main()
