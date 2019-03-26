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
    #Test default options (but with RF strandedness).
    def test_default_stranded(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        expected_annotatedjunctions = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default-stranded-annotatedjunctions.out")[0]
        expected_annotatedvariants = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default-stranded-annotatedvariants.out")[0]
        expected_junctions = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default-stranded-junctions.out")[0]
        params = ["cis-splice-effects", "identify", "-s 1",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        self.assertFilesEqual(expected_annotatedjunctions, output_annotatedjunctions, err)
        self.assertFilesEqual(expected_annotatedvariants, output_annotatedvariants, err)
        self.assertFilesEqual(expected_junctions, output_junctions, err)

    #Test default options (but with unstranded).
    def test_default(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        expected_annotatedjunctions = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default-annotatedjunctions.out")[0]
        expected_annotatedvariants = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default-annotatedvariants.out")[0]
        expected_junctions = self.inputFiles("cis-splice-effects-identify/expected-cis-splice-effects-identify-default-junctions.out")[0]
        params = ["cis-splice-effects", "identify", "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        self.assertFilesEqual(expected_annotatedjunctions, output_annotatedjunctions, err)
        self.assertFilesEqual(expected_annotatedvariants, output_annotatedvariants, err)
        self.assertFilesEqual(expected_junctions, output_junctions, err)

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

    #Test variants annotates params
    def test_exonic_intronic_distance(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        exonic_distance = "-e 6"
        intronic_distance = "-i 6"
        dont_skip_single_exon_transcripts = "-S"
        params = ["cis-splice-effects", "identify", 
                  exonic_distance, 
                  intronic_distance, 
                  dont_skip_single_exon_transcripts, 
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        #self.assertFilesEqual(expected_file, output_file, err)

    def test_allexonic(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        exonic_distance = "-E"
        intronic_distance = ""
        dont_skip_single_exon_transcripts = ""
        params = ["cis-splice-effects", "identify", 
                  exonic_distance, 
                  intronic_distance, 
                  dont_skip_single_exon_transcripts, 
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        #self.assertFilesEqual(expected_file, output_file, err)

    def test_allintronic(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        exonic_distance = ""
        intronic_distance = "-I"
        dont_skip_single_exon_transcripts = ""
        params = ["cis-splice-effects", "identify", 
                  exonic_distance, 
                  intronic_distance, 
                  dont_skip_single_exon_transcripts, 
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        #self.assertFilesEqual(expected_file, output_file, err)

    def test_allexonic_someintronic(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        exonic_distance = "-E"
        intronic_distance = "-i 6"
        dont_skip_single_exon_transcripts = ""
        params = ["cis-splice-effects", "identify", 
                  exonic_distance, 
                  intronic_distance, 
                  dont_skip_single_exon_transcripts, 
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        #self.assertFilesEqual(expected_file, output_file, err)

    def test_allintronic_someexonic(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        exonic_distance = "-e 6"
        intronic_distance = "-I"
        dont_skip_single_exon_transcripts = ""
        params = ["cis-splice-effects", "identify", 
                  exonic_distance, 
                  intronic_distance, 
                  dont_skip_single_exon_transcripts, 
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0, err)
        #self.assertFilesEqual(expected_file, output_file, err)

    #Test junctions extract params
    def test_anchor_stranded(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        print "BAM1 is ", bam1
        for anchor in ["", "30"]:
            if anchor != "":
                anchor = "-a " + anchor
            params = ["cis-splice-effects", "identify", 
                  anchor,
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
            rv, err = self.execute(params)
            self.assertEqual(rv, 0)
            #self.assertFilesEqual(expected_file, output_file)

    def test_anchor(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        print "BAM1 is ", bam1
        for anchor in ["", "30"]:
            if anchor != "":
                anchor = "-a " + anchor
            params = ["cis-splice-effects", "identify", 
                  anchor,
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
            rv, err = self.execute(params)
            self.assertEqual(rv, 0)
            #self.assertFilesEqual(expected_file, output_file)

    def test_intron_size(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_annotatedjunctions = self.tempFile("observed-cse-identify.out")
        output_annotatedvariants = self.tempFile("observed-cse-identify-variants.out")
        output_junctions = self.tempFile("observed-cse-identify-junctions.out")
        min_intron = "8039"
        max_intron = "8039"
        params = ["cis-splice-effects", "identify", 
                  "-m", min_intron, "-M", max_intron,
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)
        #self.assertFilesEqual(expected_file, output_file)

    def test_missing_bam(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = "does_not_exist.bam"
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("extract.out")
        params = ["cis-splice-effects", "identify",
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 1)

    def test_no_gtf(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        output_file = self.tempFile("extract.out")
        params = ["cis-splice-effects", "identify",
                  "-s 0",
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 1)

    def test_window_size(self):
        variants = self.inputFiles("vcf/test1.vcf")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("observed-annotate.vcf")
        variants = self.inputFiles("vcf/test1.vcf")[0]
        bam1 = self.inputFiles("bam/test_hcc1395.2.bam")[0]
        fasta = self.inputFiles("fa/test_chr22.fa")[0]
        gtf = self.inputFiles("gtf/test_ensemble_chr22.2.gtf")[0]
        output_file = self.tempFile("extract.out")
        params = ["cis-splice-effects", "identify",
                  "-s 0",
                  "-w 5"
                  "-o ", output_annotatedjunctions,
                  "-v ", output_annotatedvariants,
                  "-j ", output_junctions,
                  variants, bam1, fasta, gtf]
        rv, err = self.execute(params)
        self.assertEqual(rv, 0)

if __name__ == "__main__":
    main()
