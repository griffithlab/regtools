#!/usr/bin/env python

'''
test_junctions_extract.py -- Integration test for `regtools junctions extract`

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
'''

from integrationtest import IntegrationTest, main
import unittest

class TestExtract(IntegrationTest, unittest.TestCase):
    def test_junctions_extract_anchor_stranded(self):
        bam1 = self.inputFiles("bam/test_hcc1395.bam")[0]
        output_file = self.tempFile("extract.out")
        print "BAM1 is ", bam1
        for anchor in ["", "30"]:
            expected_file = self.inputFiles("junctions-extract/expected-stranded-a" +
                                            anchor + ".out")[0]
            if anchor != "":
                anchor = "-a " + anchor
            params = ["junctions", "extract", anchor, "-o", output_file, bam1]
            rv, err = self.execute(params)
            self.assertEqual(rv, 0)
            self.assertFilesEqual(expected_file, output_file)

    def test_junctions_extract_anchor(self):
        bam1 = self.inputFiles("bam/test_hcc1395.bam")[0]
        output_file = self.tempFile("extract.out")
        print "BAM1 is ", bam1
        for anchor in ["", "30"]:
            expected_file = self.inputFiles("junctions-extract/expected-a" +
                                            anchor + ".out")[0]
            if anchor != "":
                anchor = "-a " + anchor
            params = ["junctions", "extract", anchor, "-s 0", "-o", output_file, bam1]
            rv, err = self.execute(params)
            self.assertEqual(rv, 0)
            self.assertFilesEqual(expected_file, output_file)

    def test_junctions_extract_intron_size(self):
        bam1 = self.inputFiles("bam/test_hcc1395.bam")[0]
        output_file = self.tempFile("extract.out")
        min_intron = "8039"
        max_intron = "8039"
        expected_file = self.inputFiles("junctions-extract/expected-i" +
                min_intron + "-I" + max_intron +
                ".out")[0]
        params = ["junctions", "extract", "-s 0", "-o", output_file,
                  "-m", min_intron, "-M", max_intron, bam1]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)
        self.assertFilesEqual(expected_file, output_file)

    def test_region(self):
        bam1 = self.inputFiles("bam/test_hcc1395.bam")[0]
        output_file = self.tempFile("extract.out")
        region = "1:22405013-22405020"
        expected_file = self.inputFiles("junctions-extract/expected-r" +
                region + ".out")[0]
        params = ["junctions", "extract", "-s 0", "-o", output_file, "-r", region,
                  bam1]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)
        self.assertFilesEqual(expected_file, output_file)

    def test_no_bam(self):
        output_file = self.tempFile("extract.out")
        params = ["junctions", "extract", "-s 0", "-o", output_file]
        rv, err = self.execute(params)
        self.assertEqual(rv, 1)

    def test_help(self):
        output_file = self.tempFile("extract.out")
        params = ["junctions", "extract", "-h"]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)

if __name__ == "__main__":
    main()
