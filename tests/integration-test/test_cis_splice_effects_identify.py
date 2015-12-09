#!/usr/bin/env python

'''
test_cis_splice_effects_identify.py -- Integration test for `regtools cis-splice-effects identify`

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

class TestCisSpliceEffectsIdentify(IntegrationTest, unittest.TestCase):
    #Test default options.
    def test_default(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-cse-identify.out")
        expected_file = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default.out")[0]
        params = ["cis-splice-effects", "identify",
                  "-o ", output_file, variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        self.assertFilesEqual(expected_file, output_file, err)

    #Test -h works as expected
    def test_help(self):
        params = ["cis-splice-effects", "identify", "-h "]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)

    #Test missing input
    def test_nobam(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-cse-identify.out")
        params = ["cis-splice-effects", "identify",
                  "-o ", output_file, variants, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 1, err)

if __name__ == "__main__":
    main()
