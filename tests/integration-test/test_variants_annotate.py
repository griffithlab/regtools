#!/usr/bin/env python

'''
test_variants_annotate.py -- Integration test for `regtools variants annotate`

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

class TestAnnotate(IntegrationTest, unittest.TestCase):
    def test_variants_annotate_default(self):
        variants = self.inputFiles("test1.vcf")[0]
        gtf = self.inputFiles("test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        expected_file = self.inputFiles("variants-annotate/expected-annotate-default.out")[0]
        params = ["variants", "annotate",
                  "-o ", output_file, variants, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        self.assertFilesEqual(expected_file, output_file, err)
    def test_variants_annotate_exonic_intronic_distance(self):
        variants = self.inputFiles("test1.vcf")[0]
        gtf = self.inputFiles("test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        expected_file = self.inputFiles("variants-annotate/expected-annotate-e6-i6-S.out")[0]
        exonic_distance = "-e 6"
        intronic_distance = "-i 6"
        dont_skip_single_exon_transcripts = "-S"
        params = ["variants", "annotate", "-o ", output_file, exonic_distance,
                  intronic_distance, dont_skip_single_exon_transcripts,
                  variants, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        self.assertFilesEqual(expected_file, output_file, err)

if __name__ == "__main__":
    main()
